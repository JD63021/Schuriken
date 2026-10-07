/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_PCGNinja_H
#define Schuriken_PCGNinja_H

#include "lduMatrix.H"

namespace Foam
{

class PCGNinja
:
    public lduMatrix::solver
{
public:

    TypeName("pcgninja");

    PCGNinja
    (
        const word& fieldName,
        const lduMatrix& matrix,
        const FieldField<Field, scalar>& interfaceBouCoeffs,
        const FieldField<Field, scalar>& interfaceIntCoeffs,
        const lduInterfaceFieldPtrsList& interfaces,
        const dictionary& solverControls
    );

    virtual ~PCGNinja()
    {}

    virtual solverPerformance solve
    (
        scalarField& psi,
        const scalarField& source,
        const direction cmpt = 0
    ) const;
};

}

#endif
