/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_DistributedOps_H
#define Schuriken_DistributedOps_H

#include "label.H"
#include "scalarField.H"

#include <vector>

namespace Foam
{
namespace schuriken
{

scalar dotLocal(const scalarField& a, const scalarField& b);
scalar dotGlobal(const scalarField& a, const scalarField& b, const label comm);
scalar normGlobal(const scalarField& a, const label comm);

std::vector<scalar> dotBatchGlobal
(
    const std::vector<scalarField>& left,
    const scalarField& right,
    const label comm
);

std::vector<scalar> crossGramGlobal
(
    const std::vector<scalarField>& left,
    const std::vector<scalarField>& right,
    const label comm
);

bool normalizeEuclid
(
    scalarField& a,
    const label comm,
    const scalar tol=1e-30
);

}
}

#endif
