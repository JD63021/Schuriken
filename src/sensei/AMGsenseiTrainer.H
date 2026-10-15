/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_AMGsenseiTrainer_H
#define Schuriken_AMGsenseiTrainer_H

#include "AMGsensei.H"
#include "AMGsenseiState.H"
#include "ReducedSpace.H"
#include "lduMatrix.H"

namespace Foam
{
namespace schuriken
{

// Learn from a completed solve.  space must hold A W and W^T A W for the
// solve's matrix and current W; it is extended with W and cleared if W is
// recompressed.
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
);

void printSenseiCumulative
(
    const AMGsenseiState& state,
    const AMGsensei& config
);

}
}

#endif
