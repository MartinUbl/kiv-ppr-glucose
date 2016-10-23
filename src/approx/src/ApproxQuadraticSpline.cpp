#include "ApproxQuadraticSpline.h"

#include "../../cli/appconfig.h"
#include <amp.h>
#include <atomic>
#include <iostream>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>

// safe, no name conflict
using namespace concurrency;

void ApproxQuadraticSpline::CalculateParametersForMask(const uint32_t mask)
{
    // shifting base for mask-based calculation
    size_t base;
    // offsets of shifted values (used for fast lookup) (one extra to avoid division by 8 in loop)
    size_t offsets[9];

    size_t i;
    size_t j;

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

    j = 1;
    for (i = 1; i < valueCount - 1; i++, j++)
    {
        if (j == 8)
        {
            base += offsets[j] - offsets[0];
            j = 0;
        }

        if (base + offsets[j] > valueCount || base + offsets[j + 1] > valueCount)
            break;

        CalculateCoefsFor(mask, i, aCoefs[mask][i - 1], bCoefs[mask][i - 1],
            values[base + offsets[j]].datetime, values[base + offsets[j + 1]].datetime,
            values[base + offsets[j]].level, values[base + offsets[j + 1]].level);
    }
}

void ApproxQuadraticSpline::CalculateParameters_threads()
{
    size_t i;
    std::atomic<uint32_t> atMask = 1;

    std::thread** workers = new std::thread*[appWorkerCount];

    for (i = 0; i < appWorkerCount; i++)
    {
        workers[i] = new std::thread([&atMask](ApproxQuadraticSpline* obj) {
            uint32_t myMask;
            while ((myMask = atMask++) < APPROX_MASK_COUNT)
                obj->CalculateParametersForMask(myMask);
        }, this);
    }

    for (i = 0; i < appWorkerCount; i++)
    {
        if (workers[i]->joinable())
            workers[i]->join();

        delete workers[i];
    }

    delete workers;
}

void ApproxQuadraticSpline::CalculateParameters_AMP()
{
    // TODO
}

void ApproxQuadraticSpline::CalculateParameters_TBB()
{
    tbb::parallel_for(1, APPROX_MASK_COUNT, [&](int i) {
        CalculateParametersForMask(i);
    });
}

HRESULT IfaceCalling ApproxQuadraticSpline::Approximate(TApproximationParams *params)
{
    uint32_t mask;

    // cache count and values pointer
    mEnumeratedLevels->GetLevelsCount(&valueCount);
    mEnumeratedLevels->GetLevels(&values);

    if (appConcurrency == ConcurrencyType::ct_serial)
    {
        for (mask = 1; mask < APPROX_MASK_COUNT; mask++)
            CalculateParametersForMask(mask);
    }
    else if (appConcurrency == ConcurrencyType::ct_parallel_threads)
    {
        CalculateParameters_threads();
    }
    else if (appConcurrency == ConcurrencyType::ct_parallel_amp_gpu)
    {
        CalculateParameters_AMP();
    }
    else if (appConcurrency == ConcurrencyType::ct_parallel_tbb)
    {
        CalculateParameters_TBB();
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
    // offsets of shifted values (used for fast lookup) (one extra to avoid division by 8 in loop)
    size_t offsets[9];
    size_t i, j;

    // calculate offsets
    GetOffsetsForMask(mask, offsets);

    base = 0;

    i = 0;
    j = 0;
    floattype curtime = desiredtime;

    // fill the array
    for (; i < count; i++)
    {
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

void ApproxQuadraticSpline::CalculateCoefsFor(const uint32_t mask, size_t index, floattype aPrev, floattype bPrev, floattype xCur, floattype xNext, floattype yCur, floattype yNext)
{
    aCoefs[mask][index] = (((yCur - yNext) / (xCur - xNext)) - 2.0 * aPrev * xCur - bPrev) / (xNext - xCur);
    bCoefs[mask][index] = xCur * (2.0 * aPrev - 2.0 * aCoefs[mask][index]) + bPrev;
    cCoefs[mask][index] = yNext - aCoefs[mask][index] * xNext * xNext - bCoefs[mask][index] * xNext;
}

void ApproxQuadraticSpline::GetOffsetsForMask(const uint32_t mask, size_t* offsets)
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
