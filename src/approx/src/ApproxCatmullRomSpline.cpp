#include "ApproxCatmullRomSpline.h"

#include <iostream>
#include <thread>
#include "../../cli/appconfig.h"
#include "precalculated.h"

#include <amp.h>
#include <tbb/tbb.h>

#include "../../cli/OpenCLLoader.h"

// safe, no name conflict
using namespace concurrency;

inline floattype getVectorSquareLength(const floattype& x1, const floattype& y1, const floattype& x2, const floattype& y2)
{
    floattype dx = x2 - x1;
    floattype dy = y2 - y1;
    return dx*dx + dy*dy;
}

void ApproxCatmullRomSpline::CalculateParametersForMask(const uint32_t mask)
{
    // shifting base for mask-based calculation
    size_t base;
    // offsets of shifted values (used for fast lookup)
    const int *offsets = mask_shift_base[mask];

    const floattype tauCoef = catmullRomParam_Centripetal;

    size_t i, j;

    size_t maskedCount;
    size_t remCount;

    remCount = valueCount >> 3; //(valueCount / 8);
    maskedCount = remCount * mask_weights[mask];
    for (i = 0; i < valueCount - remCount * 8; i++)
        maskedCount += (mask >> (7 - i)) & 1;

    // resize value vector, we know how much values do we need
    xCoefs[mask].resize(maskedCount - 1);
    yCoefs[mask].resize(maskedCount - 1);

    base = 0;
    j = 0;

    size_t pos0, pos1, pos2, pos3;
    floattype dt0, dt1, dt2, t1, t2;

    for (i = 0; i < maskedCount - 2; i++, j++)
    {
        if (j == mask_weights[mask])
        {
            base += offsets[j] - offsets[0];
            j = 0;
        }

        if (base + offsets[j + 3] > valueCount)
            break;

        pos0 = base + offsets[j];
        pos1 = base + offsets[j + 1];
        pos2 = base + offsets[j + 2];
        pos3 = base + offsets[j + 3];

        dt0 = pow(getVectorSquareLength(values[pos0].datetime, values[pos0].level, values[pos1].datetime, values[pos1].level), tauCoef);
        dt1 = pow(getVectorSquareLength(values[pos1].datetime, values[pos1].level, values[pos2].datetime, values[pos2].level), tauCoef);
        dt2 = pow(getVectorSquareLength(values[pos2].datetime, values[pos2].level, values[pos3].datetime, values[pos3].level), tauCoef);

        if (dt1 < FLT_EPSILON)
            dt1 = 1.0f;
        if (dt0 < FLT_EPSILON)
            dt0 = dt1;
        if (dt2 < FLT_EPSILON)
            dt2 = dt1;

        // init xCoefs

        t1 = (values[pos1].datetime - values[pos0].datetime) / dt0 - (values[pos2].datetime - values[pos0].datetime) / (dt0 + dt1) + (values[pos2].datetime - values[pos1].datetime) / dt1;
        t2 = (values[pos2].datetime - values[pos1].datetime) / dt1 - (values[pos3].datetime - values[pos1].datetime) / (dt1 + dt2) + (values[pos3].datetime - values[pos2].datetime) / dt2;

        t1 *= dt1;
        t2 *= dt1;

        xCoefs[mask][i + 1][0] = values[pos1].datetime;
        xCoefs[mask][i + 1][1] = t1;
        xCoefs[mask][i + 1][2] = -3 * values[pos1].datetime + 3 * values[pos2].datetime - 2 * t1 - t2;
        xCoefs[mask][i + 1][3] = 2 * values[pos1].datetime - 2 * values[pos2].datetime + t1 + t2;

        // init yCoefs

        t1 = (values[pos1].level - values[pos0].level) / dt0 - (values[pos2].level - values[pos0].level) / (dt0 + dt1) + (values[pos2].level - values[pos1].level) / dt1;
        t2 = (values[pos2].level - values[pos1].level) / dt1 - (values[pos3].level - values[pos1].level) / (dt1 + dt2) + (values[pos3].level - values[pos2].level) / dt2;

        t1 *= dt1;
        t2 *= dt1;

        yCoefs[mask][i + 1][0] = values[pos1].level;
        yCoefs[mask][i + 1][1] = t1;
        yCoefs[mask][i + 1][2] = -3 * values[pos1].level + 3 * values[pos2].level - 2 * t1 - t2;
        yCoefs[mask][i + 1][3] = 2 * values[pos1].level - 2 * values[pos2].level + t1 + t2;
    }

    for (i = 0; i < 4; i++)
    {
        xCoefs[mask][0][i] = xCoefs[mask][1][i];
        yCoefs[mask][0][i] = yCoefs[mask][1][i];
    }
}

void ApproxCatmullRomSpline::CalculateParameters_threads()
{
    std::thread** workers = new std::thread*[appWorkerCount];

    // thread work assignment structure - "farmer"
    struct threadWork_t
    {
        public:
            threadWork_t(size_t* limit) { resetWork(); index_limit = limit; };

            // fills m (mask) and i (index) with current work values
            bool getWork(uint8_t &m, int32_t &i)
            {
                std::unique_lock<std::mutex> lck(mtx);

                // if we finished mask calculation
                if (index == index_limit[mask])
                {
                    // if we are at the end, no more work
                    if (mask == 255)
                        return false;

                    mask++;
                    index = 0;
                }
                else
                    index++;

                m = mask;
                i = index;
                return true;
            }

            // resets work counters
            void resetWork()
            {
                std::unique_lock<std::mutex> lck(mtx);

                mask = 1;
                index = 0;
            }
        private:
            // currently processed mask
            uint8_t mask;
            // currently processed index
            int32_t index;
            // masked index limits for masks
            size_t* index_limit;

            std::mutex mtx;
    };

    //
}

void ApproxCatmullRomSpline::CalculateParameters_AMP()
{
    //
}

void ApproxCatmullRomSpline::CalculateParameters_TBB()
{
    //
}

void ApproxCatmullRomSpline::CalculateParameters_OpenCL()
{
    //
}

HRESULT IfaceCalling ApproxCatmullRomSpline::Approximate(TApproximationParams *params)
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
        // mask-level calculation parallelization is actually FASTER than parameter-level calculation for smaller amounts of data
        // for 250 values: mask-level: ~35ms, parameter-level: ~52ms on my developer machine (i5 4570)

        // TODO: find better boundary between "small" and "big" amount of data, value of 10000 is guessed based on tendencies;
        //       such value would be most likely different for different hardware

        // for smaller amounts - mask-level calculation parallelization
        if (valueCount < 10000)
        {
            tbb::parallel_for(uint32_t(1), uint32_t(APPROX_MASK_COUNT), uint32_t(1), [&](uint32_t mask) {
                CalculateParametersForMask(mask);
            });
        }
        else // for higher amounts - parameter-level calculation parallelization
        {
            CalculateParameters_TBB();
        }
    }
    else if (appConcurrency == ConcurrencyType::ct_parallel_opencl)
    {
        CalculateParameters_OpenCL();
    }

    return S_OK;
}

HRESULT IfaceCalling ApproxCatmullRomSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder)
{
    // desired mask to be retrieved
    const uint8_t mask = appCurrentTestMask;

    HRESULT res;
    size_t index, origIndex;

    size_t remCount, maskedCount;

    remCount = (valueCount / 8);
    maskedCount = remCount * mask_weights[mask];
    for (int i = 0; i < valueCount - remCount * 8; i++)
        maskedCount += (mask >> (7 - i)) & 1;

    // retrieve index for desired time
    res = GetIndexFor(desiredtime, index, origIndex, mask);
    if (res == S_FALSE)
        return S_FALSE;

    // shifting base for mask-based calculation
    size_t base;
    size_t i, j;

    // retrieve precalculated offsets
    const int* offsets = mask_shift_base[mask];

    base = (origIndex / 8) * 8;
    if (mask_index_transform_inverse[mask][(origIndex % 8)] == 0 && base > 0)
        base -= 8;

    if (mask_index_transform_inverse[mask][origIndex % 8] == 0)
    {
        if (origIndex < 8)
            j = 0;
        else
            j = mask_index_transform_inverse[mask][7] - 1;
    }
    else
        j = mask_index_transform_inverse[mask][origIndex % 8] - 1;

    i = 0;
    floattype curtime = desiredtime;
    floattype calctime;

    // fill the array
    for (; i < count; i++)
    {
        // we want more than we have
        if (index >= maskedCount - 1)
            break;

        if (index != 0)
            calctime = (curtime - values[base + offsets[j]].datetime) / (values[base + offsets[j + 1]].datetime - values[base + offsets[j]].datetime);
        else
            calctime = (curtime - values[base + offsets[j + 1]].datetime) / (values[base + offsets[j + 1]].datetime - values[base + offsets[j]].datetime);

        // TODO: transform X coordinate (calctime) to time using xCoefs !
        // now gives slightly different results

        // no derivation: return absolute value
        if (derivationorder == 0)
            levels[i] = yCoefs[mask][index][3]*calctime*calctime*calctime + yCoefs[mask][index][2]*calctime*calctime + yCoefs[mask][index][1]*calctime + yCoefs[mask][index][0];
        // 1st order derivation
        else if (derivationorder == 1)
            levels[i] = 3*yCoefs[mask][index][3]*calctime*calctime + 2*yCoefs[mask][index][2]*calctime + yCoefs[mask][index][1];

        // move to next time
        curtime += stepping;

        // it is time to move to next index (next curve, next coefficients)
        if (curtime > values[base + offsets[j + 1]].datetime)
        {
            index++;
            // move base if needed (this is needed for mask-based calculation)
            if (index % mask_weights[mask] == 0 && index != 0)
                base += offsets[j + 1] - offsets[0];

            j = (j + 1) % mask_weights[mask];
        }
    }

    *filled = i;

    return S_OK;
}

HRESULT ApproxCatmullRomSpline::GetIndexFor(floattype time, size_t &index, size_t &origIndex, uint8_t mask)
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
            origIndex = i - 1;
            if (origIndex >= 8)
                index = (mask_index_transform_inverse[mask][7] - 1) + (mask_index_transform_inverse[mask][7] * ((origIndex / 8) - 1)) + mask_index_transform_inverse[mask][(i - 1) % 8];
            else
                index = (mask_index_transform_inverse[mask][origIndex % 8] == 0) ? 0 : (mask_index_transform_inverse[mask][origIndex % 8] - 1);

            return S_OK;
        }
    }

    return S_FALSE;
}
