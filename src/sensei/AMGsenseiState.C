/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "AMGsenseiState.H"

Foam::schuriken::AMGsenseiState::AMGsenseiState()
:
    initialized(false),
    nCells(0),
    fieldName(),
    solveCount(0),
    stockSolves(0),
    stockPCGIts(0),
    podSolves(0),
    podPCGIts(0),
    acceptedUpdates(0),
    recompressCount(0),
    W(),
    batch()
{}

void Foam::schuriken::AMGsenseiState::reset
(
    const label n,
    const word& field
)
{
    initialized = true;
    nCells = n;
    fieldName = field;
    solveCount = 0;
    stockSolves = 0;
    stockPCGIts = 0;
    podSolves = 0;
    podPCGIts = 0;
    acceptedUpdates = 0;
    recompressCount = 0;
    W.clear();
    batch.clear();
}

Foam::schuriken::AMGsenseiState& Foam::schuriken::state()
{
    static AMGsenseiState instance;
    return instance;
}
