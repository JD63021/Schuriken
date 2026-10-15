/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_DenseOps_H
#define Schuriken_DenseOps_H

#include "label.H"
#include "scalar.H"

#include <vector>

namespace Foam
{
namespace schuriken
{

std::vector<scalar> identityDense(const label n);
void symmetrize(std::vector<scalar>& A, const label n);

bool jacobiEigenSymmetric
(
    std::vector<scalar> A,
    const label n,
    std::vector<scalar>& eval,
    std::vector<scalar>& V
);

bool invertDense
(
    const std::vector<scalar>& A,
    const label n,
    std::vector<scalar>& inv
);

std::vector<scalar> denseMatVec
(
    const std::vector<scalar>& A,
    const std::vector<scalar>& x,
    const label n
);

std::vector<scalar> denseMul
(
    const std::vector<scalar>& A,
    const std::vector<scalar>& B,
    const label n
);

//- For a definite symmetric G return sign = +1 (positive definite) or -1
//  (negative definite) and Gmhalf = (sign*G)^(-1/2).  gmin/gmax are the
//  extreme eigenvalue magnitudes.  Returns false if G is indefinite or
//  numerically singular.  OpenFOAM pressure matrices are negative definite,
//  so W^T A W normally has sign = -1.
bool definiteInvSqrt
(
    const std::vector<scalar>& G,
    const label n,
    std::vector<scalar>& Gmhalf,
    scalar& sign,
    scalar& gmin,
    scalar& gmax,
    scalar& gcond
);

}
}

#endif
