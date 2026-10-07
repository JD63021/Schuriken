/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_AMGsenseiState_H
#define Schuriken_AMGsenseiState_H

#include "label.H"
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
    label solveCount;
    label stockSolves;
    label stockPCGIts;
    label podSolves;
    label podPCGIts;
    label acceptedUpdates;
    label recompressCount;
    std::vector<scalarField> W;
    std::vector<scalarField> batch;

    AMGsenseiState();

    void reset(const label n, const word& field);
};

// One in-memory adaptive state for the current pressure field.  The accessor
// keeps ownership out of the solver translation unit and gives us a clean
// place to evolve to per-region/per-field state later.
AMGsenseiState& state();

}
}

#endif
