/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_ReducedSpace_H
#define Schuriken_ReducedSpace_H

#include "label.H"
#include "scalarField.H"

#include <vector>

namespace Foam
{
namespace schuriken
{

// Products of the recycled basis W with the current matrix A: the
// distributed AW = A W and the small replicated G = W^T A W with its inverse.
struct ReducedSpace
{
    std::vector<scalarField> AW;
    std::vector<scalar> G;
    std::vector<scalar> Ginv;

    label rank() const { return static_cast<label>(AW.size()); }

    void clear()
    {
        AW.clear();
        G.clear();
        Ginv.clear();
    }
};

// Form G = W^T AW (one global reduction) and its inverse.  Returns false if
// G is singular.
bool factorReducedSpace
(
    const std::vector<scalarField>& W,
    ReducedSpace& space,
    const label comm
);

// Append v (with Av = A v) to W and border G with W^T Av and v^T Av, using
// one global reduction and no matrix products.  Returns false, leaving W
// and space unchanged, if the enlarged G is singular.  vAv returns v^T A v.
bool appendToReducedSpace
(
    std::vector<scalarField>& W,
    ReducedSpace& space,
    const scalarField& v,
    const scalarField& Av,
    scalar& vAv,
    const label comm
);

// Replace W by the numerically independent part of span(W).  Eigenvectors
// of G with non-negligible eigenvalues give A-orthogonal combinations of W;
// AW is combined the same way, so no matrix products are needed.  Clears W
// and space if nothing usable remains.  Returns the new rank.
label repairBasis
(
    std::vector<scalarField>& W,
    ReducedSpace& space,
    const label comm
);

}
}

#endif
