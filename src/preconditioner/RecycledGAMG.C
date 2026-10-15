/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "RecycledGAMG.H"
#include "AllReduce.H"
#include "DenseOps.H"
#include "DistributedOps.H"
#include "Pstream.H"

#include <cmath>

namespace Foam
{
namespace schuriken
{

void applyMatrix
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt,
    const scalarField& x,
    scalarField& Ax
)
{
    matrix.Amul(Ax, x, interfaceBouCoeffs, interfaces, cmpt);
}

bool buildReducedSpace
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt,
    const label comm,
    const std::vector<scalarField>& W,
    ReducedSpace& space
)
{
    const label r = static_cast<label>(W.size());
    const label nCells = matrix.diag().size();

    space.clear();
    space.AW.reserve(static_cast<std::size_t>(r));
    for (label j=0; j<r; ++j)
    {
        scalarField aw(nCells, 0.0);
        applyMatrix
        (
            matrix,
            interfaceBouCoeffs,
            interfaces,
            cmpt,
            W[static_cast<std::size_t>(j)],
            aw
        );
        space.AW.push_back(aw);
    }

    return factorReducedSpace(W, space, comm);
}

RecycledGAMG::RecycledGAMG
(
    lduMatrix::preconditioner& gamg,
    const std::vector<scalarField>& W,
    const ReducedSpace& space,
    const direction cmpt,
    const label comm
)
:
    gamg_(gamg),
    W_(W),
    space_(space),
    cmpt_(cmpt),
    comm_(comm),
    active_(!W.empty() && space.rank() == static_cast<label>(W.size())),
    t_(W.empty() ? 0 : W[0].size())
{}

scalar RecycledGAMG::residualSumMag
(
    const scalarField& r,
    std::vector<scalar>& Wtr
) const
{
    const label rnk = active_ ? rank() : 0;

    // [sum|r|, W^T r] in one collective
    std::vector<scalar> values(1, 0.0);
    forAll(r, cell)
    {
        values[0] += mag(r[cell]);
    }
    if (rnk)
    {
        const std::vector<scalar> local = dotBatchLocal(W_, r);
        values.insert(values.end(), local.begin(), local.end());
    }
    sumAllReduce(values.data(), rnk + 1, comm_);

    Wtr.assign(values.begin() + 1, values.end());
    return values[0];
}

scalar RecycledGAMG::apply
(
    const scalarField& r,
    const std::vector<scalar>& Wtr,
    scalarField& z
) const
{
    if (!active_)
    {
        gamg_.precondition(z, r, cmpt_);
        return gSumProd(z, r, comm_);
    }

    const label rnk = rank();

    // Coarse coefficients a = G^-1 W^T r and the residual they leave,
    // t = (I - A Q) r = r - AW a
    const std::vector<scalar> a = denseMatVec(space_.Ginv, Wtr, rnk);

    std::vector<scalar> c(a);
    for (scalar& ck : c)
    {
        ck = -ck;
    }
    t_ = r;
    addCombination(t_, space_.AW, c);

    // g = M t, written straight into z
    z = 0.0;
    gamg_.precondition(z, t_, cmpt_);

    // [(AW)^T g, r.g] in one collective
    std::vector<scalar> values = dotBatchLocal(space_.AW, z);
    values.push_back(dotLocal(r, z));
    sumAllReduce(values.data(), rnk + 1, comm_);
    const scalar rg = values.back();
    values.pop_back();

    // b = G^-1 (AW)^T g;  z = Q r + (I - Q A) g = g + W (a - b)
    // and r.z = r.g + (W^T r).(a - b) without another reduction
    const std::vector<scalar> b = denseMatVec(space_.Ginv, values, rnk);
    scalar rz = rg;
    for (label k=0; k<rnk; ++k)
    {
        const std::size_t kk = static_cast<std::size_t>(k);
        c[kk] = a[kk] - b[kk];
        rz += Wtr[kk]*c[kk];
    }
    addCombination(z, W_, c);

    return rz;
}

void RecycledGAMG::apply(const scalarField& r, scalarField& z) const
{
    std::vector<scalar> Wtr;
    if (active_)
    {
        Wtr = dotBatchGlobal(W_, r, comm_);
    }
    apply(r, Wtr, z);
}

void RecycledGAMG::gramSpectrumAbs
(
    scalar& gmin,
    scalar& gmax,
    scalar& gcond
) const
{
    gmin = GREAT;
    gmax = 0.0;
    gcond = GREAT;

    if (!active_)
    {
        return;
    }

    std::vector<scalar> eval;
    std::vector<scalar> V;
    if (jacobiEigenSymmetric(space_.G, rank(), eval, V))
    {
        for (std::size_t i=0; i<eval.size(); ++i)
        {
            gmin = min(gmin, mag(eval[i]));
            gmax = max(gmax, mag(eval[i]));
        }
        gcond = gmin > VSMALL ? gmax/gmin : GREAT;
    }
}

void RecycledGAMG::symmetryCheck
(
    scalar& relAsymP,
    scalar& relAsymM
) const
{
    const label nCells = W_.empty() ? 0 : W_[0].size();
    const label seed = Pstream::myProcNo(comm_);

    scalarField x(nCells);
    scalarField y(nCells);
    forAll(x, cell)
    {
        const scalar sx =
            std::sin(12.9898*(cell + 1) + 78.233*seed)*43758.5453;
        const scalar sy =
            std::sin(39.3468*(cell + 1) + 11.135*seed)*24634.6345;
        x[cell] = sx - std::floor(sx) - 0.5;
        y[cell] = sy - std::floor(sy) - 0.5;
    }

    const auto relAsym = [this]
    (
        const scalarField& Px,
        const scalarField& Py,
        const scalarField& x,
        const scalarField& y
    )
    {
        const scalar yPx = dotGlobal(y, Px, comm_);
        const scalar xPy = dotGlobal(x, Py, comm_);
        return mag(yPx - xPy)/max(max(mag(yPx), mag(xPy)), VSMALL);
    };

    scalarField Px(nCells, 0.0);
    scalarField Py(nCells, 0.0);
    apply(x, Px);
    apply(y, Py);
    relAsymP = relAsym(Px, Py, x, y);

    gamg_.precondition(Px, x, cmpt_);
    gamg_.precondition(Py, y, cmpt_);
    relAsymM = relAsym(Px, Py, x, y);
}

}
}
