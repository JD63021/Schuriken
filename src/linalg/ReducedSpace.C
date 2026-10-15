/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "ReducedSpace.H"
#include "DenseOps.H"
#include "DistributedOps.H"

namespace Foam
{
namespace schuriken
{

bool factorReducedSpace
(
    const std::vector<scalarField>& W,
    ReducedSpace& space,
    const label comm
)
{
    const label r = static_cast<label>(W.size());

    space.G = crossGramGlobal(W, space.AW, comm);
    symmetrize(space.G, r);
    if (!invertDense(space.G, r, space.Ginv))
    {
        return false;
    }

    // A symmetric inverse keeps the balanced preconditioner symmetric
    symmetrize(space.Ginv, r);
    return true;
}

bool appendToReducedSpace
(
    std::vector<scalarField>& W,
    ReducedSpace& space,
    const scalarField& v,
    const scalarField& Av,
    scalar& vAv,
    const label comm
)
{
    const label r = space.rank();
    const label r1 = r + 1;

    // [W^T Av; v^T Av] in one reduction
    W.push_back(v);
    const std::vector<scalar> g = dotBatchGlobal(W, Av, comm);

    std::vector<scalar> G(static_cast<std::size_t>(r1*r1), 0.0);
    for (label i=0; i<r; ++i)
    {
        for (label j=0; j<r; ++j)
        {
            G[static_cast<std::size_t>(i*r1+j)] =
                space.G[static_cast<std::size_t>(i*r+j)];
        }
        G[static_cast<std::size_t>(i*r1+r)] = g[static_cast<std::size_t>(i)];
        G[static_cast<std::size_t>(r*r1+i)] = g[static_cast<std::size_t>(i)];
    }
    G[static_cast<std::size_t>(r*r1+r)] = g[static_cast<std::size_t>(r)];

    std::vector<scalar> Ginv;
    if (!invertDense(G, r1, Ginv))
    {
        W.pop_back();
        return false;
    }
    symmetrize(Ginv, r1);

    space.AW.push_back(Av);
    space.G.swap(G);
    space.Ginv.swap(Ginv);
    vAv = g[static_cast<std::size_t>(r)];
    return true;
}

label repairBasis
(
    std::vector<scalarField>& W,
    ReducedSpace& space,
    const label comm
)
{
    const label r = static_cast<label>(W.size());

    std::vector<scalar> eval;
    std::vector<scalar> V;
    if (r == 0 || !jacobiEigenSymmetric(space.G, r, eval, V))
    {
        W.clear();
        space.clear();
        return 0;
    }

    // The dominant eigenvalue fixes the sign of a definite A
    scalar emax = 0.0;
    scalar sign = 1.0;
    for (label k=0; k<r; ++k)
    {
        const scalar e = eval[static_cast<std::size_t>(k)];
        if (mag(e) > emax)
        {
            emax = mag(e);
            sign = e > 0.0 ? 1.0 : -1.0;
        }
    }

    const label nCells = W[0].size();
    std::vector<scalarField> Wnew;
    std::vector<scalarField> AWnew;

    for (label k=0; k<r; ++k)
    {
        if (!(sign*eval[static_cast<std::size_t>(k)] > 1e-10*emax))
        {
            continue;
        }

        scalarField w(nCells, 0.0);
        scalarField aw(nCells, 0.0);
        for (label j=0; j<r; ++j)
        {
            const scalar vjk = V[static_cast<std::size_t>(j*r+k)];
            const scalarField& wj = W[static_cast<std::size_t>(j)];
            const scalarField& awj = space.AW[static_cast<std::size_t>(j)];
            forAll(w, cell)
            {
                w[cell] += vjk*wj[cell];
                aw[cell] += vjk*awj[cell];
            }
        }

        const scalar nw = normGlobal(w, comm);
        if (!(nw > VSMALL))
        {
            continue;
        }
        w /= nw;
        aw /= nw;
        Wnew.push_back(w);
        AWnew.push_back(aw);
    }

    W.swap(Wnew);
    space.AW.swap(AWnew);

    if (W.empty() || !factorReducedSpace(W, space, comm))
    {
        W.clear();
        space.clear();
    }

    return static_cast<label>(W.size());
}

}
}
