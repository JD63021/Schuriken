/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_AllReduce_H
#define Schuriken_AllReduce_H

#include "label.H"
#include "scalar.H"

namespace Foam
{
namespace schuriken
{

// In-place global sum of n values over communicator comm.  Built against
// MPI this is a single MPI_Allreduce; OpenFOAM-12 itself only exposes
// scalar all-reduces, otherwise a list gather/scatter is used.
void sumAllReduce(scalar* values, const label n, const label comm);

}
}

#endif
