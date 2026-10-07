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

    const auto appendAOrthogonal =
    [&](scalarField v) -> bool
    {
        if (!normalizeEuclid(v, comm))
        {
            return false;
        }

        for (label pass=0; pass<2; ++pass)
        {
            if (state.W.empty())
            {
                break;
            }

            std::vector<scalarField> curAW;
            std::vector<scalar> curG;
            std::vector<scalar> curGinv;
            if
            (
                !buildReducedSpace
                (
                    matrix,
                    interfaceBouCoeffs,
                    interfaces,
                    cmpt,
                    comm,
                    state.W,
                    curAW,
                    curG,
                    curGinv
                )
            )
            {
                return false;
            }

            scalarField Av(nCells, 0.0);
            applyMatrix
            (
                matrix,
                interfaceBouCoeffs,
                interfaces,
                cmpt,
                v,
                Av
            );

            const label r = static_cast<label>(state.W.size());
            const std::vector<scalar> h = dotBatchGlobal(state.W, Av, comm);
            const std::vector<scalar> c = denseMatVec(curGinv, h, r);

            for (label j=0; j<r; ++j)
            {
                const scalar cj = c[static_cast<std::size_t>(j)];
                const scalarField& wj = state.W[static_cast<std::size_t>(j)];
                forAll(v, cell)
                {
                    v[cell] -= cj*wj[cell];
                }
            }
        }

        const scalar projectedNorm = normGlobal(v, comm);
        if (!normalizeEuclid(v, comm) || !(projectedNorm > 1e-10))
        {
            return false;
        }

        scalarField Av(nCells, 0.0);
        applyMatrix
        (
            matrix,
            interfaceBouCoeffs,
            interfaces,
            cmpt,
            v,
            Av
        );
        const scalar aEnergy = dotGlobal(v, Av, comm);

        state.W.push_back(v);
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

    std::vector<scalarField> fullAW;
    std::vector<scalar> fullG;
    std::vector<scalar> fullGinv;
    const label oldR = static_cast<label>(state.W.size());

    if
    (
        !buildReducedSpace
        (
            matrix,
            interfaceBouCoeffs,
            interfaces,
            cmpt,
            comm,
            state.W,
            fullAW,
            fullG,
            fullGinv
        )
    )
    {
        return;
    }

    std::vector<scalar> Gmhalf;
    scalar gmin = 0.0;
    scalar gmax = 0.0;
    scalar gcond = 0.0;

    if (!spdInvSqrt(fullG, oldR, Gmhalf, gmin, gmax, gcond))
    {
        Info
            << "SCHURIKEN_AMGSENSEI_RECOMPRESS_SKIP"
            << " solve=" << solveNo
            << " rank=" << oldR
            << " reason=G_NOT_SPD"
            << " status=PASS" << endl;
        return;
    }

    std::vector<scalarField> BAW;
    BAW.reserve(static_cast<std::size_t>(oldR));
    for (label j=0; j<oldR; ++j)
    {
        scalarField baw(nCells, 0.0);
        gamg.precondition(baw, fullAW[static_cast<std::size_t>(j)], cmpt);
        BAW.push_back(baw);
    }

    std::vector<scalar> K = crossGramGlobal(fullAW, BAW, comm);
    symmetrize(K, oldR);

    std::vector<scalar> R = denseMul(denseMul(Gmhalf, K, oldR), Gmhalf, oldR);
    symmetrize(R, oldR);

    std::vector<scalar> mu;
    std::vector<scalar> Y;
    if (!jacobiEigenSymmetric(R, oldR, mu, Y))
    {
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

    const label keep = min(recompressRank, oldR);
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
        for (label j=0; j<oldR; ++j)
        {
            const scalar cj = coeff[static_cast<std::size_t>(j)];
            const scalarField& wj = state.W[static_cast<std::size_t>(j)];
            forAll(wnew, cell)
            {
                wnew[cell] += cj*wj[cell];
            }
        }

        if (!normalizeEuclid(wnew, comm))
        {
            FatalErrorInFunction
                << "Spectral recompression produced zero/nonfinite mode"
                << exit(FatalError);
        }
        Wnew.push_back(wnew);
    }

    state.W.swap(Wnew);
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
        << " solves=" << state.solveCount
        << " rank=" << state.W.size()
        << " stockSolves=" << state.stockSolves
        << " stockPCGIts=" << state.stockPCGIts
        << " podSolves=" << state.podSolves
        << " podPCGIts=" << state.podPCGIts
        << " totalPCGIts=" << (state.stockPCGIts + state.podPCGIts)
        << " acceptedUpdates=" << state.acceptedUpdates
        << " recompressCount=" << state.recompressCount
        << " learnOddOnly=" << config.learnOddOnly()
        << " enabled=" << (config.enabled() ? 1 : 0)
        << " status=PASS" << endl;
}

}
}
