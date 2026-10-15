/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "PCGNinja.H"
#include "AMGsensei.H"
#include "AMGsenseiState.H"
#include "AMGsenseiTrainer.H"
#include "RecycledGAMG.H"
#include "ReducedSpace.H"
#include "Pstream.H"

#include <chrono>
#include <vector>

namespace Foam
{
    defineTypeNameAndDebug(PCGNinja, 0);

    lduMatrix::solver::
    addsymMatrixConstructorToTable<PCGNinja>
        addPCGNinjaSymMatrixConstructorToTable_;
}

namespace
{
    // Monotonic seconds; the wall clock may jump (e.g. under WSL)
    double monotonicSeconds()
    {
        return std::chrono::duration<double>
        (
            std::chrono::steady_clock::now().time_since_epoch()
        ).count();
    }
}

Foam::PCGNinja::PCGNinja
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}

Foam::solverPerformance Foam::PCGNinja::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    using schuriken::AMGsenseiState;
    using schuriken::RecycledGAMG;

    const label nCells = psi.size();
    const label comm = matrix().mesh().comm();

    const AMGsensei sensei(controlDict_);
    const bool senseiEnabled = sensei.enabled();
    const label trainSnapshots = sensei.trainSnapshots();
    const label addRank = sensei.addRank();
    const label maxRank = sensei.maxRank();
    const label recompressRank = sensei.recompressRank();
    const label tailIter = sensei.tailIter();
    const label learnOddOnly = sensei.learnOddOnly();

    if
    (
        trainSnapshots < 1
     || addRank < 1
     || maxRank < addRank
     || recompressRank < 1
     || recompressRank >= maxRank
     || tailIter < 1
    )
    {
        FatalErrorInFunction
            << "Bad AMGsensei controls: trainSnapshots=" << trainSnapshots
            << " addRank=" << addRank
            << " maxRank=" << maxRank
            << " recompressRank=" << recompressRank
            << " tailIter=" << tailIter
            << exit(FatalError);
    }

    AMGsenseiState& senseiState =
        schuriken::state(matrix().mesh(), fieldName_, sensei.persist());

    if (senseiState.restoredRank >= 0 && senseiState.solveCount == 0)
    {
        Info
            << "SCHURIKEN_AMGSENSEI_RESTART"
            << " region=" << senseiState.region
            << " field=" << fieldName_
            << " rank=" << senseiState.restoredRank
            << " cellsFile=" << senseiState.nCells
            << " cellsMesh=" << nCells
            << " status="
            << (senseiState.nCells == nCells ? "PASS" : "REJECTED")
            << endl;
    }

    if
    (
        !senseiState.initialized
     || senseiState.nCells != nCells
     || senseiState.fieldName != fieldName_
    )
    {
        senseiState.reset(nCells, fieldName_);
        Info
            << "SCHURIKEN_AMGSENSEI_RESET"
            << " region=" << senseiState.region
            << " field=" << fieldName_
            << " cells=" << nCells
            << " status=PASS" << endl;
    }

    ++senseiState.solveCount;
    const label solveNo = senseiState.solveCount;
    const bool learnThisSolve = sensei.learnThisSolve(solveNo);

    if (solveNo == 1)
    {
        Info
            << "SCHURIKEN_AMGSENSEI_CONFIG"
            << " region=" << senseiState.region
            << " field=" << fieldName_
            << " enabled=" << (senseiEnabled ? 1 : 0)
            << " trainSnapshots=" << trainSnapshots
            << " addRank=" << addRank
            << " maxRank=" << maxRank
            << " recompressRank=" << recompressRank
            << " recompress=" << (sensei.recompress() ? "yes" : "no")
            << " persist=" << (sensei.persist() ? "yes" : "no")
            << " tailIter=" << tailIter
            << " learnOddOnly=" << learnOddOnly
            << " learn=" << sensei.learnMode()
            << " parallel=" << (Pstream::parRun() ? 1 : 0)
            << " nProcs=" << Pstream::nProcs()
            << " reductionPolicy=BATCHED_GLOBAL_POD"
            << " status=PASS" << endl;
    }

    double tick = monotonicSeconds();
    const auto lap = [&tick]()
    {
        const double now = monotonicSeconds();
        const double dt = now - tick;
        tick = now;
        return dt;
    };

    autoPtr<lduMatrix::preconditioner> gamgPtr =
        lduMatrix::preconditioner::New(*this, controlDict_);

    // A disabled AMGsensei is a pure stock PCG+GAMG path, even if an earlier
    // solve populated the in-memory basis before runtime dictionary changes.
    const std::vector<scalarField> emptyBasis;
    const std::vector<scalarField>& activeBasis =
        senseiEnabled ? senseiState.W : emptyBasis;

    // A singular W^T A W must not stop the CFD run: keep the independent
    // part of the basis, or fall back to stock PCG+GAMG with an empty one.
    schuriken::ReducedSpace space;
    if
    (
        !activeBasis.empty()
     && !schuriken::buildReducedSpace
        (
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            cmpt,
            comm,
            activeBasis,
            space
        )
    )
    {
        const label rankOld = senseiState.W.size();
        const label rankNew =
            schuriken::repairBasis(senseiState.W, space, comm);

        Info
            << "SCHURIKEN_AMGSENSEI_REPAIR"
            << " solve=" << solveNo
            << " rankOld=" << rankOld
            << " rankNew=" << rankNew
            << " reason=G_SINGULAR"
            << " status=PASS" << endl;
    }

    RecycledGAMG recycled(*gamgPtr, activeBasis, space, cmpt, comm);

    const bool usedPOD = recycled.active();
    if (usedPOD)
    {
        scalar gmin = GREAT;
        scalar gmax = 0.0;
        scalar gcond = GREAT;
        recycled.gramSpectrumAbs(gmin, gmax, gcond);

        Info
            << "SCHURIKEN_AMGSENSEI_SETUP"
            << " solve=" << solveNo
            << " rank=" << recycled.rank()
            << " GminAbs=" << gmin
            << " GmaxAbs=" << gmax
            << " GcondAbs=" << gcond
            << " policy=CURRENT_A"
            << " status=PASS" << endl;

        if (sensei.checkSymmetry())
        {
            scalar relAsymP = 0;
            scalar relAsymM = 0;
            recycled.symmetryCheck(relAsymP, relAsymM);
            Info
                << "SCHURIKEN_AMGSENSEI_SYMMETRY"
                << " solve=" << solveNo
                << " rank=" << recycled.rank()
                << " relAsymP=" << relAsymP
                << " relAsymM=" << relAsymM
                << " status=PASS" << endl;
        }
    }

    const double setupDt = lap();
    senseiState.setupSeconds += setupDt;

    solverPerformance solverPerf(word("pcgninja"), fieldName_);

    const scalarField psiStart(psi);
    scalarField psiTail(nCells, 0.0);
    bool tailSaved = false;

    scalarField pA(nCells);
    scalarField wA(nCells);
    scalar wArA = solverPerf.great_;
    scalar wArAold = wArA;

    schuriken::applyMatrix
    (
        matrix_,
        interfaceBouCoeffs_,
        interfaces_,
        cmpt,
        psi,
        wA
    );
    scalarField rA(source - wA);

    const scalar normFactor = this->normFactor(psi, source, wA, pA);

    // Every iteration needs three collectives, as stock PCG: W^T r travels
    // with the residual norm and (AW)^T M t with r.z
    std::vector<scalar> Wtr;
    solverPerf.initialResidual() =
        recycled.residualSumMag(rA, Wtr)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    if
    (
        minIter_ > 0
     || !solverPerf.checkConvergence(tolerance_, relTol_)
    )
    {
        do
        {
            wArAold = wArA;
            wArA = recycled.apply(rA, Wtr, wA);

            if (solverPerf.nIterations() == 0)
            {
                pA = wA;
            }
            else
            {
                const scalar beta = wArA/wArAold;
                forAll(pA, cell)
                {
                    pA[cell] = wA[cell] + beta*pA[cell];
                }
            }

            schuriken::applyMatrix
            (
                matrix_,
                interfaceBouCoeffs_,
                interfaces_,
                cmpt,
                pA,
                wA
            );

            const scalar wApA = gSumProd(wA, pA, matrix().mesh().comm());
            if (solverPerf.checkSingularity(mag(wApA)/normFactor))
            {
                break;
            }

            const scalar alpha = wArA/wApA;
            forAll(psi, cell)
            {
                psi[cell] += alpha*pA[cell];
                rA[cell] -= alpha*wA[cell];
            }

            const label completedIter = solverPerf.nIterations() + 1;
            if (learnThisSolve && completedIter == tailIter)
            {
                psiTail = psi;
                tailSaved = true;
                Info
                    << "SCHURIKEN_AMGSENSEI_TAIL_ITERATE"
                    << " solve=" << solveNo
                    << " pcgIter=" << completedIter
                    << " source=X_AT_FIXED_KRYLOV_DEPTH"
                    << " status=PASS" << endl;
            }

            solverPerf.finalResidual() =
                recycled.residualSumMag(rA, Wtr)/normFactor;
        }
        while
        (
            (
                ++solverPerf.nIterations() < maxIter_
             && !solverPerf.checkConvergence(tolerance_, relTol_)
            )
         || solverPerf.nIterations() < minIter_
        );
    }

    const label finalIts = solverPerf.nIterations();
    const double iterDt = lap();
    senseiState.iterSeconds += iterDt;

    if (usedPOD)
    {
        ++senseiState.podSolves;
        senseiState.podPCGIts += finalIts;
    }
    else
    {
        ++senseiState.stockSolves;
        senseiState.stockPCGIts += finalIts;
    }

    Info
        << "SCHURIKEN_AMGSENSEI_SOLVE"
        << " region=" << senseiState.region
        << " solve=" << solveNo
        << " learn=" << (learnThisSolve ? 1 : 0)
        << " rank=" << senseiState.W.size()
        << " usedPOD=" << (usedPOD ? 1 : 0)
        << " pCG=" << finalIts
        << " initialResidual=" << solverPerf.initialResidual()
        << " finalResidual=" << solverPerf.finalResidual()
        << " seconds=" << setupDt + iterDt
        << " status=PASS" << endl;

    schuriken::updateSenseiAfterSolve
    (
        senseiState,
        sensei,
        solveNo,
        learnThisSolve,
        psi,
        psiStart,
        psiTail,
        tailSaved,
        finalIts,
        matrix_,
        interfaceBouCoeffs_,
        interfaces_,
        *gamgPtr,
        space,
        cmpt,
        comm
    );
    senseiState.trainSeconds += lap();

    schuriken::printSenseiCumulative(senseiState, sensei);

    return solverPerf;
}
