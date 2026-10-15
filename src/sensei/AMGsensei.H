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

    // Explicit "learn" switch: 1 yes, 0 no, -1 unset (use learnOddOnly)
    label learn_;

    // Spectral recompression at maxRank; otherwise the full basis is kept
    // and learning stops
    bool recompress_;

    // Write the basis at write times and read it back on restart
    bool persist_;

    bool checkSymmetry_;

public:

    explicit AMGsensei(const dictionary& solverControls);

    bool enabled() const { return enabled_; }
    label trainSnapshots() const { return trainSnapshots_; }
    label addRank() const { return addRank_; }
    label maxRank() const { return maxRank_; }
    label recompressRank() const { return recompressRank_; }
    bool recompress() const { return recompress_; }
    bool persist() const { return persist_; }
    label tailIter() const { return tailIter_; }
    label learnOddOnly() const { return learnOddOnly_; }

    //- "yes", "no" or "auto" (learn switch not given)
    word learnMode() const
    {
        return learn_ < 0 ? "auto" : (learn_ ? "yes" : "no");
    }

    //- Whether the solveNo-th solve of the field trains the basis: the
    //  explicit learn switch of this solver dictionary if given, otherwise
    //  every solve, or every odd solve with learnOddOnly
    bool learnThisSolve(const label solveNo) const
    {
        if (!enabled_)
        {
            return false;
        }
        if (learn_ >= 0)
        {
            return learn_;
        }
        return !learnOddOnly_ || solveNo % 2 == 1;
    }

    //- Test diagnostic: compare y'Px with x'Py in every recycled solve
    bool checkSymmetry() const { return checkSymmetry_; }
};

}

#endif
