#include "CommonApproxSupport.h"
#include "precalculated.h"
#include <tbb\tbb.h>

floattype Cubic_IdentitySolve(floattype a, floattype b, floattype c, floattype d)
{
    floattype re[3], im[3];

    if (a == 0 || d == 0)
        return 0;

    b /= a;
    c /= a;
    d /= a;
    floattype disc, q, r, dum1, s, t, term1, r13;

    q = (3.0*c - (b*b)) / 9.0;
    r = -(27.0*d / 54.0) + b*(9.0*c / 54.0 - 2.0*(b*b) / 54.0);

    disc = q*q*q + r*r;
    im[0] = 0;
    term1 = (b / 3.0);
    // only one real solution, two complex
    if (disc > 0)
    {
        s = r + sqrt(disc);
        s = ((s < 0) ? -pow(-s, (1.0 / 3.0)) : pow(s, (1.0 / 3.0)));
        t = r - sqrt(disc);
        t = ((t < 0) ? -pow(-t, (1.0 / 3.0)) : pow(t, (1.0 / 3.0)));
        re[0] = -term1 + s + t;

        // if the solution is between 0 and 1, accept it
        if (re[0] >= 0.0 && re[0] <= 1.0)
            return re[0];
    }
    else // all three real solutions
    {
        im[2] = im[1] = 0;
        if (disc == 0)
        {
            r13 = ((r < 0) ? -pow(-r, (1.0 / 3.0)) : pow(r, (1.0 / 3.0)));
            re[0] = -term1 + 2.0*r13;
            re[2] = re[1] = -(r13 + term1);
        }
        else
        {
            q = -q;
            dum1 = acos(r / sqrt(q*q*q));
            r13 = 2.0*sqrt(q);
            re[0] = -term1 + r13*cos(dum1 / 3.0);

            re[1] = -term1 + r13*cos((dum1 + 2.0*M_PI) / 3.0);
            re[2] = -term1 + r13*cos((dum1 + 4.0*M_PI) / 3.0);
        }

        for (int i = 0; i < 3; i++)
        {
            if (re[i] > 0.0 && re[i] <= 1.0)
                return re[i];
        }
    }

    return 0;
}

CFindMaskedBounds::CFindMaskedBounds(TGlucoseLevel* levels, uint8_t mask, uint32_t skipSides[2]) : mLevels(levels), mMask(mask), mSkipSides{ skipSides[0], skipSides[1] }
{
    mBounds.MaxLevel = mLevels[mask_shift_base[mMask][0 + skipSides[0]]].level;
    mBounds.MaxTime = mLevels[mask_shift_base[mMask][0 + skipSides[0]]].datetime;
    mBounds.MinLevel = mBounds.MaxLevel;
    mBounds.MinTime = mBounds.MaxTime;
}

CFindMaskedBounds::CFindMaskedBounds(CFindMaskedBounds& x, tbb::split) : mLevels(x.mLevels), mBounds(x.mBounds), mMask(x.mMask), mSkipSides{x.mSkipSides[0], mSkipSides[1]}
{
    //
}

void CFindMaskedBounds::operator()(const tbb::blocked_range<size_t>& r)
{
    TGlucoseLevel *levels = mLevels;
    TGlucoseLevelBounds bounds = mBounds;
    size_t end = r.end();
    size_t pos;

    for (size_t i = r.begin(); i != end; ++i)
    {
        if (i < mSkipSides[0] || i >= mSkipSides[1])
            continue;

        pos = 8 * (i / mask_weights[mMask]) + mask_index_transform[mMask][i % mask_weights[mMask]];

        bounds.MaxLevel = max(bounds.MaxLevel, levels[pos].level);
        bounds.MinLevel = min(bounds.MinLevel, levels[pos].level);

        bounds.MaxTime = max(bounds.MaxTime, levels[pos].datetime);
        bounds.MinTime = min(bounds.MinTime, levels[pos].datetime);
    }

    mBounds = bounds;
}

void CFindMaskedBounds::join(const CFindMaskedBounds& y)
{
    mBounds.MaxLevel = max(mBounds.MaxLevel, y.mBounds.MaxLevel);
    mBounds.MinLevel = min(mBounds.MinLevel, y.mBounds.MinLevel);

    mBounds.MaxTime = max(mBounds.MaxTime, y.mBounds.MaxTime);
    mBounds.MinTime = min(mBounds.MinTime, y.mBounds.MinTime);
}
