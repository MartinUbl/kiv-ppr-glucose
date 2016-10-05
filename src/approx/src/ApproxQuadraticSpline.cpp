#include "ApproxQuadraticSpline.h"

#include <iostream>

HRESULT IfaceCalling ApproxQuadraticSpline::Approximate(TApproximationParams *params)
{
    // shifting base for mask-based calculation
    size_t base;
    // offsets of shifted values (used for fast lookup)
    size_t offsets[8];

    size_t i;
    uint32_t mask;

    // cache count and values pointer
    mEnumeratedLevels->GetLevelsCount(&valueCount);
    mEnumeratedLevels->GetLevels(&values);

    // calculate parameters for all mask
    for (mask = 0; mask < APPROX_MASK_COUNT; mask++)
    {
        // resize value vectors, we know how much values do we need
        aCoefs[mask].resize(valueCount - 1);
        bCoefs[mask].resize(valueCount - 1);
        cCoefs[mask].resize(valueCount - 1);

        // calculate offsets
        GetOffsetsForMask(mask, offsets);

        base = 0;

        // initial condition for quadratic spline
        aCoefs[mask][0] = 0;
        bCoefs[mask][0] = (values[offsets[1]].level - values[offsets[0]].level) / (values[offsets[1]].datetime - values[offsets[0]].datetime);
        cCoefs[mask][0] = values[offsets[0]].level - bCoefs[mask][0] * values[offsets[0]].datetime;

        for (i = 1; i < valueCount - 1; i++)
        {
            if (base + offsets[i] > valueCount || base + offsets[i + 1] > valueCount)
                break;

            CalculateCoefsFor(mask, i, aCoefs[mask][i - 1], bCoefs[mask][i - 1], values[base + offsets[i]].datetime, values[base + offsets[i + 1]].datetime, values[base + offsets[i]].level, values[base + offsets[i + 1]].level);

            if (i % 8 == 0)
                base += 8;
        }
    }

    return S_OK;
}

HRESULT IfaceCalling ApproxQuadraticSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder)
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
    // offsets of shifted values (used for fast lookup)
    size_t offsets[8];
    size_t i;

    // calculate offsets
    GetOffsetsForMask(mask, offsets);

    base = 0;

    i = 0;
    floattype curtime = desiredtime;

    // fill the array
    for (; i < count; i++)
    {
        // it is time to move to next index (next curve, next coefficients)
        if (index < 7 && curtime > values[base + offsets[index + 1]].datetime)
            index++;

        // we want more than we have
        if (index >= valueCount)
            break;

        // no derivation: return absolute value
        if (derivationorder == 0)
            levels[i] = curtime*(aCoefs[mask][index]*curtime + bCoefs[mask][index]) + cCoefs[mask][index];
        // 1st order derivation, quadratic function derives as 2*a*x + b
        else if (derivationorder == 1)
            levels[i] = 2 * aCoefs[mask][index] * curtime + bCoefs[mask][index];

        // move to next time
        curtime += stepping;

        // move base if needed (this is needed for mask-based calculation)
        if (index % 8 == 0 && index != 0)
            base += 8;
    }

    *filled = i - 1;

    return S_OK;
}

void ApproxQuadraticSpline::CalculateCoefsFor(const uint32_t mask, size_t index, floattype aPrev, floattype bPrev, floattype xCur, floattype xNext, floattype yCur, floattype yNext)
{
    aCoefs[mask][index] = (((yCur - yNext) / (xCur - xNext)) - 2.0 * aPrev * xCur - bPrev) / (xNext - xCur);
    bCoefs[mask][index] = xCur * (2.0 * aPrev - 2.0 * aCoefs[mask][index]) + bPrev;
    cCoefs[mask][index] = yNext - aCoefs[mask][index] * xNext * xNext - bCoefs[mask][index] * xNext;
}

void ApproxQuadraticSpline::GetOffsetsForMask(const uint32_t mask, size_t* offsets)
{
    size_t i, j;

    // initial state is: 0, 1, 2, 3, 4, 5, 6, 7 - this applies to full mask
    for (i = 0; i < 8; i++)
        offsets[i] = i;

    for (i = 0; i < 8; i++)
    {
        // for each 0 bit in mask we shift current and following offsets by one
        if (((1 << i) & mask) == 0)
        {
            for (j = i; j < 8; j++)
                offsets[j]++;
        }
    }
}

HRESULT ApproxQuadraticSpline::GetIndexFor(floattype time, size_t &index)
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
