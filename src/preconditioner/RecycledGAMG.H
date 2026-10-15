/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_RecycledGAMG_H
#define Schuriken_RecycledGAMG_H

#include "lduMatrix.H"
#include "ReducedSpace.H"

#include <vector>

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
);

// Form AW = A W and factor G = W^T A W.  Returns false if G is singular.
bool buildReducedSpace
(
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt,
    const label comm,
    const std::vector<scalarField>& W,
    ReducedSpace& space
);

class RecycledGAMG
{
    lduMatrix::preconditioner& gamg_;
    const std::vector<scalarField>& W_;
    const ReducedSpace& space_;
    direction cmpt_;
    label comm_;
    bool active_;

    //- Work field for apply()
    mutable scalarField t_;

public:

    // W and space must describe the same basis for the current matrix
    RecycledGAMG
    (
        lduMatrix::preconditioner& gamg,
        const std::vector<scalarField>& W,
        const ReducedSpace& space,
        const direction cmpt,
        const label comm
    );

    bool active() const { return active_; }
    label rank() const { return static_cast<label>(W_.size()); }

    //- Global sum |r| for the residual norm and, when recycling, the global
    //  W^T r needed by the next apply, reduced together in one collective
    scalar residualSumMag(const scalarField& r, std::vector<scalar>& Wtr) const;

    //- z = P r given Wtr = W^T r from residualSumMag.  Returns the global
    //  r.z, which is reduced together with (AW)^T M t in one collective
    scalar apply
    (
        const scalarField& r,
        const std::vector<scalar>& Wtr,
        scalarField& z
    ) const;

    //- z = P r, forming W^T r itself (diagnostics)
    void apply(const scalarField& r, scalarField& z) const;

    void gramSpectrumAbs
    (
        scalar& gmin,
        scalar& gmax,
        scalar& gcond
    ) const;

    //- Relative asymmetry |y'Px - x'Py|/max(|y'Px|, |x'Py|) of the recycled
    //  preconditioner P and of the base preconditioner M for fixed
    //  pseudo-random x and y
    void symmetryCheck(scalar& relAsymP, scalar& relAsymM) const;
};

}
}

#endif
