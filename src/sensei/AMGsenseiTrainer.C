/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "AMGsenseiTrainer.H"
#include "DenseOps.H"
#include "DistributedOps.H"
#include "PODBasis.H"
#include "RecycledGAMG.H"

#include <algorithm>
#include <vector>

namespace Foam
{
namespace schuriken
{

void updateSenseiAfterSolve
(
    AMGsenseiState& state,
    const AMGsensei& config,
    const label solveNo,
    const bool learnThisSolve,
    const scalarField& psi,
    const scalarField& psiStart,
    const scalarField& psiTail,
    const bool tailSaved,
    const label finalIts,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    lduMatrix::preconditioner& gamg,
    ReducedSpace& space,
    const direction cmpt,
    const label comm
)
{
    if (!learnThisSolve)
    {
        return;
    }

    const label nCells = psi.size();
    const label trainSnapshots = config.trainSnapshots();
    const label addRank = config.addRank();
    const label maxRank = config.maxRank();
    const label recompressRank = config.recompressRank();
    const label tailIter = config.tailIter();

    // Without recompression a full basis is kept as it is
    if
    (
        !config.recompress()
     && static_cast<label>(state.W.size()) >= maxRank
    )
    {
        state.batch.clear();
        return;
    }

    scalarField snapshot(nCells, 0.0);
    bool haveSnapshot = false;
    const bool seedMode = state.W.empty();

    if (seedMode)
    {
        forAll(snapshot, cell)
        {
            snapshot[cell] = psi[cell] - psiStart[cell];
        }
        haveSnapshot = normalizeEuclid(snapshot, comm);
        if (haveSnapshot)
        {
            Info
                << "SCHURIKEN_AMGSENSEI_SNAPSHOT"
                << " solve=" << solveNo
                << " snapshot=" << (state.batch.size()+1)
                << "/" << trainSnapshots
                << " kind=SEED_CORRECTION"
                << " status=PASS" << endl;
        }
    }
    else if (tailSaved && finalIts > tailIter)
    {
        forAll(snapshot, cell)
        {
            snapshot[cell] = psi[cell] - psiTail[cell];
        }
        const scalar rawNorm = normGlobal(snapshot, comm);
        haveSnapshot = normalizeEuclid(snapshot, comm);
        if (haveSnapshot)
        {
            Info
                << "SCHURIKEN_AMGSENSEI_TAIL_CORRECTION"
                << " solve=" << solveNo
                << " pcgIterSaved=" << tailIter
                << " finalIts=" << finalIts
                << " remainingIts=" << (finalIts-tailIter)
                << " snapshot=" << (state.batch.size()+1)
                << "/" << trainSnapshots
                << " rawNorm=" << rawNorm
                << " source=X_FINAL_MINUS_X_TAILITER"
                << " status=PASS" << endl;
        }
    }

    if (haveSnapshot)
    {
        state.batch.push_back(snapshot);
    }

    if (static_cast<label>(state.batch.size()) < trainSnapshots)
    {
        return;
    }

    std::vector<scalarField> modes;
    scalar retainedEnergy = 0.0;
    if (!snapshotPOD(state.batch, addRank, modes, retainedEnergy, comm))
    {
        state.batch.clear();
        return;
    }

    state.batch.clear();

    // The solve's reduced space normally describes W for this matrix; rebuild
    // it only if it does not (e.g. learning without a recycled solve).
    if
    (
        space.rank() != static_cast<label>(state.W.size())
     && !buildReducedSpace
        (
            matrix,
            interfaceBouCoeffs,
            interfaces,
            cmpt,
            comm,
            state.W,
            space
        )
    )
    {
        const label rankOld = state.W.size();
        repairBasis(state.W, space, comm);
        Info
            << "SCHURIKEN_AMGSENSEI_REPAIR"
            << " solve=" << solveNo
            << " rankOld=" << rankOld
            << " rankNew=" << state.W.size()
            << " reason=G_SINGULAR_TRAINER"
            << " status=PASS" << endl;
    }

    // W and its reduced space grow together.  The matrix is unchanged since
    // the solve, so only the new column needs matrix products.
    const auto appendAOrthogonal =
    [&](scalarField v) -> bool
    {
        if (!normalizeEuclid(v, comm))
        {
            return false;
        }

        scalarField Av(nCells, 0.0);
        applyMatrix(matrix, interfaceBouCoeffs, interfaces, cmpt, v, Av);

        // Two passes of A-orthogonal projection; AW updates A v alongside v
        const label r = space.rank();
        for (label pass=0; pass<2 && r>0; ++pass)
        {
            std::vector<scalar> c =
                denseMatVec(space.Ginv, dotBatchGlobal(state.W, Av, comm), r);
            for (scalar& cj : c)
            {
                cj = -cj;
            }
            addCombination(v, state.W, c);
            addCombination(Av, space.AW, c);
        }

        const scalar projectedNorm = normGlobal(v, comm);
        if (!normalizeEuclid(v, comm) || !(projectedNorm > 1e-10))
        {
            return false;
        }

        // Exact A v for the stored column
        applyMatrix(matrix, interfaceBouCoeffs, interfaces, cmpt, v, Av);

        scalar aEnergy = 0.0;
        if (!appendToReducedSpace(state.W, space, v, Av, aEnergy, comm))
        {
            return false;
        }

        Info
            << "SCHURIKEN_AMGSENSEI_CANDIDATE"
            << " solve=" << solveNo
            << " rankNow=" << state.W.size()
            << " projectedEuclidNorm=" << projectedNorm
            << " Aenergy=" << aEnergy
            << " orthogonalization=CURRENT_A"
            << " status=PASS" << endl;
        return true;
    };

    const label oldRank = static_cast<label>(state.W.size());
    label added = 0;
    for (std::size_t i=0; i<modes.size(); ++i)
    {
        if (static_cast<label>(state.W.size()) >= maxRank)
        {
            break;
        }
        if (appendAOrthogonal(modes[i]))
        {
            ++added;
        }
    }

    if (added > 0)
    {
        ++state.acceptedUpdates;
        Info
            << "SCHURIKEN_AMGSENSEI_AUGMENT"
            << " solve=" << solveNo
            << " rankOld=" << oldRank
            << " rankNew=" << state.W.size()
            << " add=" << added
            << " batchRetainedEnergy=" << retainedEnergy
            << " noProbe=1"
            << " status=PASS" << endl;
    }

    if (static_cast<label>(state.W.size()) < maxRank)
    {
        return;
    }

    if (!config.recompress())
    {
        Info
            << "SCHURIKEN_AMGSENSEI_BASIS_FULL"
            << " solve=" << solveNo
            << " rank=" << state.W.size()
            << " recompress=no learning=stopped"
            << " status=PASS" << endl;
        return;
    }

    const label oldR = static_cast<label>(state.W.size());
    const label keep = min(recompressRank, oldR);

    // If the spectral selection is impossible, keep the newest vectors so
    // the basis can still adapt instead of freezing at maxRank.
    const auto keepNewest = [&](const char* reason)
    {
        std::vector<scalarField> Wnew
        (
            state.W.end() - keep,
            state.W.end()
        );
        state.W.swap(Wnew);
        space.clear();
        ++state.recompressCount;

        Info
            << "SCHURIKEN_AMGSENSEI_RECOMPRESS"
            << " solve=" << solveNo
            << " cycle=" << state.recompressCount
            << " rankOld=" << oldR
            << " rankNew=" << state.W.size()
            << " method=keep_newest"
            << " reason=" << reason
            << " status=PASS" << endl;
    };

    if (space.rank() != oldR)
    {
        keepNewest("SPACE_MISMATCH");
        return;
    }

    std::vector<scalar> Gmhalf;
    scalar gSign = 1.0;
    scalar gmin = 0.0;
    scalar gmax = 0.0;
    scalar gcond = 0.0;

    if (!definiteInvSqrt(space.G, oldR, Gmhalf, gSign, gmin, gmax, gcond))
    {
        keepNewest("G_NOT_DEFINITE");
        return;
    }

    std::vector<scalarField> BAW;
    BAW.reserve(static_cast<std::size_t>(oldR));
    for (label j=0; j<oldR; ++j)
    {
        scalarField baw(nCells, 0.0);
        gamg.precondition(baw, space.AW[static_cast<std::size_t>(j)], cmpt);
        BAW.push_back(baw);
    }

    // K = (AW)^T M (AW) has the same sign as G for a definite A and M, so
    // scaling both by gSign keeps the generalized Ritz values mu positive.
    std::vector<scalar> K = crossGramGlobal(space.AW, BAW, comm);
    symmetrize(K, oldR);
    for (scalar& k : K)
    {
        k *= gSign;
    }

    std::vector<scalar> R = denseMul(denseMul(Gmhalf, K, oldR), Gmhalf, oldR);
    symmetrize(R, oldR);

    std::vector<scalar> mu;
    std::vector<scalar> Y;
    if (!jacobiEigenSymmetric(R, oldR, mu, Y))
    {
        keepNewest("RITZ_EIGEN_FAILED");
        return;
    }

    std::vector<label> order(static_cast<std::size_t>(oldR));
    for (label i=0; i<oldR; ++i)
    {
        order[static_cast<std::size_t>(i)] = i;
    }
    std::sort
    (
        order.begin(),
        order.end(),
        [&](const label a, const label b)
        {
            return mu[static_cast<std::size_t>(a)]
                 < mu[static_cast<std::size_t>(b)];
        }
    );

    const scalar muMin = mu[static_cast<std::size_t>(order.front())];
    const scalar muMax = mu[static_cast<std::size_t>(order.back())];
    const scalar muSelMax =
        mu[static_cast<std::size_t>(order[static_cast<std::size_t>(keep-1)])];
    const scalar muDropMin =
        keep < oldR
      ? mu[static_cast<std::size_t>(order[static_cast<std::size_t>(keep)])]
      : muSelMax;

    std::vector<scalarField> Wnew;
    Wnew.reserve(static_cast<std::size_t>(keep));

    for (label kk=0; kk<keep; ++kk)
    {
        const label ec = order[static_cast<std::size_t>(kk)];
        std::vector<scalar> coeff(static_cast<std::size_t>(oldR), 0.0);

        for (label j=0; j<oldR; ++j)
        {
            scalar s = 0.0;
            for (label a=0; a<oldR; ++a)
            {
                s +=
                    Gmhalf[static_cast<std::size_t>(j*oldR+a)]
                   *Y[static_cast<std::size_t>(a*oldR+ec)];
            }
            coeff[static_cast<std::size_t>(j)] = s;
        }

        scalarField wnew(nCells, 0.0);
        addCombination(wnew, state.W, coeff);

        if (!normalizeEuclid(wnew, comm))
        {
            keepNewest("ZERO_OR_NONFINITE_MODE");
            return;
        }
        Wnew.push_back(wnew);
    }

    // space no longer describes W; the next solve rebuilds it for its matrix
    state.W.swap(Wnew);
    space.clear();
    ++state.recompressCount;

    Info
        << "SCHURIKEN_AMGSENSEI_AMG_SPECTRUM"
        << " solve=" << solveNo
        << " rankOld=" << oldR
        << " select=" << keep
        << " muMin=" << muMin
        << " muSelectedMax=" << muSelMax
        << " muDroppedMin=" << muDropMin
        << " muMax=" << muMax
        << " spectralCondition="
        << (mag(muMin) > VSMALL ? muMax/mag(muMin) : GREAT)
        << " Gcond=" << gcond
        << " Gsign=" << gSign
        << " metric=GENERALIZED_RITZ_AWB0AW_OVER_WAW"
        << " status=PASS" << endl;

    Info
        << "SCHURIKEN_AMGSENSEI_RECOMPRESS"
        << " solve=" << solveNo
        << " cycle=" << state.recompressCount
        << " rankOld=" << oldR
        << " rankNew=" << state.W.size()
        << " method=amg_spectral"
        << " noProbe=1"
        << " status=PASS" << endl;
}

void printSenseiCumulative
(
    const AMGsenseiState& state,
    const AMGsensei& config
)
{
    Info
        << "SCHURIKEN_AMGSENSEI_CUMULATIVE"
        << " region=" << state.region
        << " field=" << state.fieldName
        << " solves=" << state.solveCount
        << " rank=" << state.W.size()
        << " stockSolves=" << state.stockSolves
        << " stockPCGIts=" << state.stockPCGIts
        << " podSolves=" << state.podSolves
        << " podPCGIts=" << state.podPCGIts
        << " totalPCGIts=" << (state.stockPCGIts + state.podPCGIts)
        << " acceptedUpdates=" << state.acceptedUpdates
        << " recompressCount=" << state.recompressCount
        << " setupSeconds=" << state.setupSeconds
        << " iterSeconds=" << state.iterSeconds
        << " trainSeconds=" << state.trainSeconds
        << " learnOddOnly=" << config.learnOddOnly()
        << " enabled=" << (config.enabled() ? 1 : 0)
        << " status=PASS" << endl;
}

}
}
