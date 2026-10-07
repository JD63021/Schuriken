/* Schuriken - OpenFOAM Foundation v12 pressure-solver add-on
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#ifndef Schuriken_AMGsensei_H
#define Schuriken_AMGsensei_H

#include "dictionary.H"
#include "label.H"

namespace Foam
{

// Small configuration wrapper for the adaptive POD trainer/recycler used by
// pcgninja.  The numerical state remains solver-internal ; this class
// deliberately owns only the public configuration surface.
class AMGsensei
{
    bool enabled_;
    label trainSnapshots_;
    label addRank_;
    label maxRank_;
    label recompressRank_;
    label tailIter_;
    label learnOddOnly_;

public:

    explicit AMGsensei(const dictionary& solverControls);

    bool enabled() const { return enabled_; }
    label trainSnapshots() const { return trainSnapshots_; }
    label addRank() const { return addRank_; }
    label maxRank() const { return maxRank_; }
    label recompressRank() const { return recompressRank_; }
    label tailIter() const { return tailIter_; }
    label learnOddOnly() const { return learnOddOnly_; }
};

}

#endif
