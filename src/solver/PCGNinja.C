/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "PCGNinja.H"
#include "AMGsensei.H"
#include "AMGsenseiState.H"
#include "AMGsenseiTrainer.H"
#include "RecycledGAMG.H"
#include "Pstream.H"

#include <vector>

namespace Foam
{
    defineTypeNameAndDebug(PCGNinja, 0);

    lduMatrix::solver::
    addsymMatrixConstructorToTable<PCGNinja>
        addPCGNinjaSymMatrixConstructorToTable_;
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

    AMGsenseiState& senseiState = schuriken::state();

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
            << " field=" << fieldName_
            << " cells=" << nCells
            << " status=PASS" << endl;
    }

    ++senseiState.solveCount;
    const label solveNo = senseiState.solveCount;
    const bool learnThisSolve =
        senseiEnabled && (!learnOddOnly || (solveNo % 2 == 1));

    if (solveNo == 1)
    {
        Info
            << "SCHURIKEN_AMGSENSEI_CONFIG"
            << " enabled=" << (senseiEnabled ? 1 : 0)
            << " trainSnapshots=" << trainSnapshots
            << " addRank=" << addRank
            << " maxRank=" << maxRank
            << " recompressRank=" << recompressRank
            << " tailIter=" << tailIter
            << " learnOddOnly=" << learnOddOnly
            << " parallel=" << (Pstream::parRun() ? 1 : 0)
            << " nProcs=" << Pstream::nProcs()
            << " reductionPolicy=BATCHED_GLOBAL_POD"
            << " status=PASS" << endl;
    }

    autoPtr<lduMatrix::preconditioner> gamgPtr =
        lduMatrix::preconditioner::New(*this, controlDict_);

    // A disabled AMGsensei is a pure stock PCG+GAMG path, even if an earlier
    // solve populated the in-memory basis before runtime dictionary changes.
    const std::vector<scalarField> emptyBasis;
    const std::vector<scalarField>& activeBasis =
        senseiEnabled ? senseiState.W : emptyBasis;

    RecycledGAMG recycled
    (
        matrix_,
        interfaceBouCoeffs_,
        interfaces_,
        *gamgPtr,
        activeBasis,
        cmpt,
        comm
    );

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
    }

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

    solverPerf.initialResidual() =
        gSumMag(rA, matrix().mesh().comm())/normFactor;
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
            recycled.apply(rA, wA);

            wArA = gSumProd(wA, rA, matrix().mesh().comm());

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
                gSumMag(rA, matrix().mesh().comm())/normFactor;
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
        << " solve=" << solveNo
        << " learn=" << (learnThisSolve ? 1 : 0)
        << " rank=" << senseiState.W.size()
        << " usedPOD=" << (usedPOD ? 1 : 0)
        << " pCG=" << finalIts
        << " initialResidual=" << solverPerf.initialResidual()
        << " finalResidual=" << solverPerf.finalResidual()
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
        cmpt,
        comm
    );

    schuriken::printSenseiCumulative(senseiState, sensei);

    return solverPerf;
}
