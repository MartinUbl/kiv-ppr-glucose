#include "ApproxAkimaSpline.h"

#include <iostream>
#include <thread>

HRESULT IfaceCalling ApproxAkimaSpline::Approximate(TApproximationParams *params)
{
    // shifting base for mask-based calculation
    size_t base;
    // offsets of shifted values (used for fast lookup) (one extra to avoid division by 8 in loop)
    size_t offsets[9];

    size_t i;
    size_t j;
    uint32_t mask;

    // cache count and values pointer
    mEnumeratedLevels->GetLevelsCount(&valueCount);
    mEnumeratedLevels->GetLevels(&values);

    std::vector<floattype> derivations;

    size_t maskedCount;
    size_t remCount;

    floattype tmpWeight1, tmpWeight2;
    floattype* m;

    // calculate parameters for all mask
    for (mask = 1; mask < APPROX_MASK_COUNT; mask++)
    {
        maskedCount = 0;
        for (i = 0; i < 8; i++)
            maskedCount += (mask >> i) & 1;

        remCount = (valueCount / 8);

        maskedCount = remCount * maskedCount;

        for (i = 0; i < valueCount - remCount * 8; i++)
            maskedCount += (mask >> (7 - i)) & 1;

        // resize value vectors, we know how much values do we need
        a0Coefs[mask].resize(valueCount - 1);
        a1Coefs[mask].resize(valueCount - 1);
        a2Coefs[mask].resize(valueCount - 1);
        a3Coefs[mask].resize(valueCount - 1);

        derivations.clear();
        derivations.resize(valueCount + 4);

        m = &derivations[2];

        // calculate offsets
        GetOffsetsForMask(mask, offsets);

        base = 0;
        j = 0;

        // precalculate derivatives
        for (i = 0; i < valueCount - 1; i++, j++)
        {
            if (j == 8)
            {
                base += offsets[j] - offsets[0];
                j = 0;
            }

            if (base + offsets[j] > valueCount || base + offsets[j + 1] > valueCount)
                break;

            m[i] = (values[base + offsets[j + 1]].level - values[base + offsets[j]].level) /
                (values[base + offsets[j + 1]].datetime - values[base + offsets[j]].datetime);
        }
        m[-1] = 2 * m[0] - m[1];
        m[-2] = 2 * m[-1] - m[0];
        m[maskedCount - 1] = 2 * m[maskedCount - 2] - m[maskedCount - 3];
        m[maskedCount] = 2 * m[maskedCount - 1] - m[maskedCount - 2];

        base = 0;
        j = 0;

        for (i = 0; i < valueCount - 1; i++, j++)
        {
            if (j == 8)
            {
                base += offsets[j] - offsets[0];
                j = 0;
            }

            if (base + offsets[j] > valueCount || base + offsets[j + 1] > valueCount)
                break;

            const floattype wsum1 = fabs(m[i + 1] - m[i]) + fabs(m[i - 1] - m[i - 2]);
            if (wsum1 == 0)
            {
                a1Coefs[mask][i] = m[i];
                a2Coefs[mask][i] = 0;
                a3Coefs[mask][i] = 0;
            }
            else
            {
                const floattype deltax = values[base + offsets[j + 1]].datetime - values[base + offsets[j]].datetime;

                // TODO: weighting coefficients!

                a0Coefs[mask][i] = values[base + offsets[j]].level;
                a1Coefs[mask][i] = m[i];
                a2Coefs[mask][i] = (3 * m[i] - 2 * a1Coefs[mask][i] - m[i + 1]) / deltax;
                a3Coefs[mask][i] = (a1Coefs[mask][i] + m[i + 1] - 2 * m[i]) / (deltax * deltax);
            }
        }
    }

    return S_OK;
}

HRESULT IfaceCalling ApproxAkimaSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder)
{
    // desired mask to be retrieved, full mask for now
    // this should not be here after testing stage
    const uint8_t mask = 0xFF;

    HRESULT res;
    size_t index;

    // retrieve index for desired time
    res = GetIndexFor(desiredtime, index);
    if (res == S_FALSE)
        return S_FALSE;

    // shifting base for mask-based calculation
    size_t base;
    // offsets of shifted values (used for fast lookup) (one extra to avoid division by 8 in loop)
    size_t offsets[9];
    size_t i, j;

    // calculate offsets
    GetOffsetsForMask(mask, offsets);

    base = 0;

    i = 0;
    j = 0;
    floattype curtime = desiredtime;
    floattype timediff;

    // fill the array
    for (; i < count; i++)
    {
        // we want more than we have
        if (index >= valueCount)
            break;

        timediff = curtime - values[base + offsets[j]].datetime;

        // no derivation: return absolute value
        if (derivationorder == 0)
            levels[i] = a0Coefs[mask][index] + timediff * (a1Coefs[mask][index] + timediff * (a2Coefs[mask][index] + timediff * a3Coefs[mask][index]));
        // 1st order derivation
        else if (derivationorder == 1)
            levels[i] = a1Coefs[mask][index] + timediff * (2.0 * a2Coefs[mask][index] + 3 * a3Coefs[mask][index] * timediff);

        // move to next time
        curtime += stepping;

        // it is time to move to next index (next curve, next coefficients)
        if (curtime > values[base + offsets[j + 1]].datetime)
        {
            index++;
            // move base if needed (this is needed for mask-based calculation)
            if (index % 8 == 0 && index != 0)
                base += offsets[j + 1] - offsets[0];

            j = (j + 1) % 8;
        }
    }

    *filled = i - 1;

    return S_OK;
}

void ApproxAkimaSpline::GetOffsetsForMask(const uint32_t mask, size_t* offsets)
{
    size_t i, j;

    i = 0;
    j = 0;
    // calculate one extra offset to avoid division in main loop
    while (i < 9)
    {
        // inverse the mask, so it matches the group in "the right direction" (MSB = 0. value, LSB = 7. value within group)
        if (((1 << (7 - (j % 8))) & mask) != 0)
        {
            offsets[i] = j;
            i++;
        }

        j++;
    }
}

HRESULT ApproxAkimaSpline::GetIndexFor(floattype time, size_t &index)
{
    // secure out of range values
    if (time + DBL_EPSILON < values[0].datetime || time > values[valueCount - 1].datetime)
        return S_FALSE;

    // for now, linear search is implemented as there is higher chance searched value would be at the beginning
    // since this method is called only at the beginning of curve sampling

    size_t i = 1;
    for (; i < valueCount; i++)
    {
        if (time < values[i].datetime)
        {
            // this is safe since we guarantee the bounds with previous conditions
            index = i - 1;
            return S_OK;
        }
    }

    return S_FALSE;
}
