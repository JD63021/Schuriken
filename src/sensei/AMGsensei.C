/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "AMGsensei.H"

Foam::AMGsensei::AMGsensei(const dictionary& solverControls)
:
    enabled_(true),
    trainSnapshots_(3),
    addRank_(2),
    maxRank_(24),
    recompressRank_(16),
    tailIter_(10),
    learnOddOnly_(1)
{
    // Preferred Schuriken API.
    if (solverControls.found("AMGsensei"))
    {
        const dictionary& d = solverControls.subDict("AMGsensei");

        enabled_ = d.lookupOrDefault<bool>("enabled", enabled_);
        trainSnapshots_ =
            d.lookupOrDefault<label>("trainSnapshots", trainSnapshots_);
        addRank_ = d.lookupOrDefault<label>("addRank", addRank_);
        maxRank_ = d.lookupOrDefault<label>("maxRank", maxRank_);
        recompressRank_ =
            d.lookupOrDefault<label>("recompressRank", recompressRank_);
        tailIter_ = d.lookupOrDefault<label>("tailIter", tailIter_);
        learnOddOnly_ =
            d.lookupOrDefault<label>("learnOddOnly", learnOddOnly_);
    }
    else
    {
        // Compatibility with the research prototype dictionaries.
        trainSnapshots_ = solverControls.lookupOrDefault<label>
        (
            "podTrainSnapshots", trainSnapshots_
        );
        addRank_ = solverControls.lookupOrDefault<label>
        (
            "podAddRank", addRank_
        );
        maxRank_ = solverControls.lookupOrDefault<label>
        (
            "podMaxRank", maxRank_
        );
        recompressRank_ = solverControls.lookupOrDefault<label>
        (
            "podRecompressRank", recompressRank_
        );
        tailIter_ = solverControls.lookupOrDefault<label>
        (
            "podTailIter", tailIter_
        );
        learnOddOnly_ = solverControls.lookupOrDefault<label>
        (
            "podLearnOddOnly", learnOddOnly_
        );
    }
}
