/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_AMGsenseiState_H
#define Schuriken_AMGsenseiState_H

#include "label.H"
#include "lduMesh.H"
#include "scalarField.H"
#include "word.H"

#include <vector>

namespace Foam
{
namespace schuriken
{

struct AMGsenseiState
{
    bool initialized;
    label nCells;
    word fieldName;

    // Mesh region the state belongs to ("none" outside an object registry)
    word region;

    label solveCount;
    label stockSolves;
    label stockPCGIts;
    label podSolves;
    label podPCGIts;
    label acceptedUpdates;
    label recompressCount;

    // Monotonic wall time spent in pcgninja: basis setup, PCG iterations
    // and training
    scalar setupSeconds;
    scalar iterSeconds;
    scalar trainSeconds;

    // Rank of a basis read back on restart, -1 if none was read
    label restoredRank;

    std::vector<scalarField> W;
    std::vector<scalarField> batch;

    AMGsenseiState();

    void reset(const label n, const word& field);
};

// Adaptive state of one field on one mesh.  It is held in the mesh's object
// registry, so every region and field learns separately and the state is
// destroyed with its mesh.  With persist the basis is written as
// <time>/AMGsensei_<field> at write times and read back from the start time
// when the state is first created.  A matrix whose mesh is not a registry
// (e.g. a GAMG coarse level) uses a process-wide table keyed by field name.
AMGsenseiState& state
(
    const lduMesh& mesh,
    const word& fieldName,
    const bool persist
);

}
}

#endif
