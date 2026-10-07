/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_RecycledGAMG_H
#define Schuriken_RecycledGAMG_H

#include "lduMatrix.H"

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
);

class RecycledGAMG
{
    const lduMatrix& matrix_;
    const FieldField<Field, scalar>& interfaceBouCoeffs_;
    const lduInterfaceFieldPtrsList& interfaces_;
    lduMatrix::preconditioner& gamg_;
    const std::vector<scalarField>& W_;
    direction cmpt_;
    label comm_;

    std::vector<scalarField> AW_;
    std::vector<scalar> G_;
    std::vector<scalar> Ginv_;
    bool active_;

public:

    RecycledGAMG
    (
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        lduMatrix::preconditioner& gamg,
        const std::vector<scalarField>& W,
        const direction cmpt,
        const label comm
    );

    bool active() const { return active_; }
    label rank() const { return static_cast<label>(W_.size()); }

    const std::vector<scalarField>& AW() const { return AW_; }
    const std::vector<scalar>& G() const { return G_; }
    const std::vector<scalar>& Ginv() const { return Ginv_; }

    void apply(const scalarField& r, scalarField& z) const;

    void gramSpectrumAbs
    (
        scalar& gmin,
        scalar& gmax,
        scalar& gcond
    ) const;
};

}
}

#endif
