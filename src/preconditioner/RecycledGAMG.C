/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "RecycledGAMG.H"
#include "DenseOps.H"
#include "DistributedOps.H"

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
    std::vector<scalarField>& AW,
    std::vector<scalar>& G,
    std::vector<scalar>& Ginv
)
{
    const label r = static_cast<label>(W.size());
    const label nCells = matrix.diag().size();

    AW.clear();
    AW.reserve(static_cast<std::size_t>(r));
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
        AW.push_back(aw);
    }

    G = crossGramGlobal(W, AW, comm);
    symmetrize(G, r);
    return invertDense(G, r, Ginv);
}

RecycledGAMG::RecycledGAMG
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    lduMatrix::preconditioner& gamg,
    const std::vector<scalarField>& W,
    const direction cmpt,
    const label comm
)
:
    matrix_(matrix),
    interfaceBouCoeffs_(interfaceBouCoeffs),
    interfaces_(interfaces),
    gamg_(gamg),
    W_(W),
    cmpt_(cmpt),
    comm_(comm),
    AW_(),
    G_(),
    Ginv_(),
    active_(!W.empty())
{
    if
    (
        active_
     && !buildReducedSpace
        (
            matrix_,
            interfaceBouCoeffs_,
            interfaces_,
            cmpt_,
            comm_,
            W_,
            AW_,
            G_,
            Ginv_
        )
    )
    {
        FatalErrorInFunction
            << "Current reduced matrix W^T A W is singular. rank="
            << W_.size()
            << exit(FatalError);
    }
}

void RecycledGAMG::apply(const scalarField& r, scalarField& z) const
{
    if (!active_)
    {
        gamg_.precondition(z, r, cmpt_);
        return;
    }

    const label rnk = rank();
    const label nCells = r.size();

    const std::vector<scalar> rhs = dotBatchGlobal(W_, r, comm_);
    const std::vector<scalar> a = denseMatVec(Ginv_, rhs, rnk);

    scalarField q(nCells, 0.0);
    scalarField t(r);
    for (label k=0; k<rnk; ++k)
    {
        const scalar ak = a[static_cast<std::size_t>(k)];
        const scalarField& wk = W_[static_cast<std::size_t>(k)];
        const scalarField& awk = AW_[static_cast<std::size_t>(k)];
        forAll(q, cell)
        {
            q[cell] += ak*wk[cell];
            t[cell] -= ak*awk[cell];
        }
    }

    scalarField g(nCells, 0.0);
    gamg_.precondition(g, t, cmpt_);

    const std::vector<scalar> rhs2 = dotBatchGlobal(AW_, g, comm_);
    const std::vector<scalar> b = denseMatVec(Ginv_, rhs2, rnk);

    z = q;
    forAll(z, cell)
    {
        z[cell] += g[cell];
    }

    for (label k=0; k<rnk; ++k)
    {
        const scalar bk = b[static_cast<std::size_t>(k)];
        const scalarField& wk = W_[static_cast<std::size_t>(k)];
        forAll(z, cell)
        {
            z[cell] -= bk*wk[cell];
        }
    }
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
    if (jacobiEigenSymmetric(G_, rank(), eval, V))
    {
        for (std::size_t i=0; i<eval.size(); ++i)
        {
            gmin = min(gmin, mag(eval[i]));
            gmax = max(gmax, mag(eval[i]));
        }
        gcond = gmin > VSMALL ? gmax/gmin : GREAT;
    }
}

}
}
