/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_PODBasis_H
#define Schuriken_PODBasis_H

#include "label.H"
#include "scalarField.H"

#include <vector>

namespace Foam
{
namespace schuriken
{

bool snapshotPOD
(
    const std::vector<scalarField>& snapshots,
    const label requestedRank,
    std::vector<scalarField>& modes,
    scalar& retainedEnergy,
    const label comm
);

}
}

#endif
