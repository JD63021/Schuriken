/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "PODBasis.H"
#include "DenseOps.H"
#include "DistributedOps.H"

#include <algorithm>
#include <cmath>

namespace Foam
{
namespace schuriken
{

bool snapshotPOD
(
    const std::vector<scalarField>& snapshots,
    const label requestedRank,
    std::vector<scalarField>& modes,
    scalar& retainedEnergy,
    const label comm
)
{
    const label m = static_cast<label>(snapshots.size());
    if (m <= 0)
    {
        retainedEnergy = 0.0;
        modes.clear();
        return false;
    }

    const label n = snapshots[0].size();
    std::vector<scalar> C = crossGramGlobal(snapshots, snapshots, comm);
    symmetrize(C, m);

    std::vector<scalar> eval;
    std::vector<scalar> V;
    if (!jacobiEigenSymmetric(C, m, eval, V))
    {
        return false;
    }

    std::vector<label> order(static_cast<std::size_t>(m));
    for (label i=0; i<m; ++i) order[static_cast<std::size_t>(i)] = i;
    std::sort
    (
        order.begin(), order.end(),
        [&](const label a, const label b)
        {
            return eval[static_cast<std::size_t>(a)] > eval[static_cast<std::size_t>(b)];
        }
    );

    scalar total = 0.0;
    for (label i=0; i<m; ++i)
    {
        if (eval[static_cast<std::size_t>(i)] > 0.0)
        {
            total += eval[static_cast<std::size_t>(i)];
        }
    }

    const scalar lam0 = Foam::max(eval[static_cast<std::size_t>(order[0])], scalar(0));
    scalar kept = 0.0;
    modes.clear();

    const label want = Foam::min(requestedRank, m);
    for (label kk=0; kk<want; ++kk)
    {
        const label ec = order[static_cast<std::size_t>(kk)];
        const scalar lam = eval[static_cast<std::size_t>(ec)];
        if (!(lam > 1e-12*Foam::max(lam0, scalar(1e-300))))
        {
            continue;
        }

        scalarField u(n, 0.0);
        const scalar invSig = 1.0/std::sqrt(lam);
        for (label q=0; q<m; ++q)
        {
            const scalar a = V[static_cast<std::size_t>(q*m+ec)]*invSig;
            const scalarField& sq = snapshots[static_cast<std::size_t>(q)];
            forAll(u, i)
            {
                u[i] += a*sq[i];
            }
        }
        if (normalizeEuclid(u, comm))
        {
            modes.push_back(u);
            kept += lam;
        }
    }

    retainedEnergy = total > 0.0 ? kept/total : 0.0;
    return !modes.empty();
}

}
}
