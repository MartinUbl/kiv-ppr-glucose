#pragma once

#include "../../common/iface/ApproxIface.h"
#include <tbb/tbb.h>

#ifndef M_PI
// Ludolph's number - PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Solve cubic equation with only one root in range <0;1>
 */
floattype Cubic_IdentitySolve(floattype a, floattype b, floattype c, floattype d);

/*
 * Class used for TBB to find bounds of masked values
 */
class CFindMaskedBounds
{
    public:
        // default "parent" constructor
        CFindMaskedBounds(TGlucoseLevel* levels, uint8_t mask, uint32_t skipSides[2]);

        // split constructor for TBB
        CFindMaskedBounds(CFindMaskedBounds& x, tbb::split);

        // function operator to make functor - find bounds in range
        void operator()(const tbb::blocked_range<size_t>& r);

        // join function to merge results
        void join(const CFindMaskedBounds& y);

        // calculated bounds
        TGlucoseLevelBounds mBounds;

    private:
        // stored levels
        TGlucoseLevel *mLevels;
        // mask used for calculation
        uint8_t mMask;
        // left and right limitation of index assigned for current run
        uint32_t mSkipSides[2];
};
