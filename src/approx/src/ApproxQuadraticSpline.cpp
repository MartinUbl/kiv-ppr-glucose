#include "ApproxQuadraticSpline.h"
#include "CommonApproxSupport.h"

#include "../../cli/appconfig.h"
#include "precalculated.h"
#include <amp.h>
#include <atomic>
#include <iostream>
#include <tbb/tbb.h>
#include <tbb/parallel_for.h>
#include "../../cli/OpenCLLoader.h"

// safe, no name conflict
using namespace concurrency;

void ApproxQuadraticSpline::CalculateParametersForMask(const uint32_t mask)
{
    // shifting base for mask-based calculation
    size_t base;

    size_t i;
    size_t j;

    // resize value vectors, we know how much values do we need
    aCoefs[mask].resize(valueCount - 1);
    bCoefs[mask].resize(valueCount - 1);
    cCoefs[mask].resize(valueCount - 1);

    // retrieve offsets
    const int* offsets = mask_shift_base[mask];

    base = 0;

    // initial condition for quadratic spline
    aCoefs[mask][0] = 0;
    bCoefs[mask][0] = (values[offsets[1]].level - values[offsets[0]].level) / (values[offsets[1]].datetime - values[offsets[0]].datetime);
    cCoefs[mask][0] = values[offsets[0]].level - bCoefs[mask][0] * values[offsets[0]].datetime;

    j = 1;
    for (i = 1; i < valueCount - 1; i++, j++)
    {
        if (j == mask_weights[mask])
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

    delete[] workers;
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

void ApproxQuadraticSpline::CalculateParameters_OpenCL()
{
    cl_int ret;

    clProgramRecord* program = GetCLProgramRecord(apxmQuadraticSpline);
    if (!program)
        return;

    void* vals_cp;
    int vals_size;

    // we have to copy values to intermediate array since we are not sure the GPU supports double precision
    // also we cannot "copy" double to float at binary level due to different sizes and schemas
    // and finally, we use this monstrosity to avoid passing array of structure instances to OpenCL program
    if (clSupportsDouble())
    {
        vals_size = sizeof(double);
        vals_cp = new double[valueCount * 2];
        double* vals_cp_typed = reinterpret_cast<double*>(vals_cp);
        for (int i = 0; i < valueCount; i++)
        {
            vals_cp_typed[i * 2 + 0] = (double)values[i].datetime;
            vals_cp_typed[i * 2 + 1] = (double)values[i].level;
        }
    }
    else
    {
        vals_size = sizeof(float);
        vals_cp = new float[valueCount * 2];
        float* vals_cp_typed = reinterpret_cast<float*>(vals_cp);
        for (int i = 0; i < valueCount; i++)
        {
            vals_cp_typed[i * 2 + 0] = (float)values[i].datetime;
            vals_cp_typed[i * 2 + 1] = (float)values[i].level;
        }
    }

    // read-only parameters
    cl_mem mask_weights_m = clCreateBuffer(program->context, CL_MEM_READ_ONLY, APPROX_MASK_COUNT * sizeof(int), NULL, &ret);
    cl_mem mask_shift_base_m = clCreateBuffer(program->context, CL_MEM_READ_ONLY, APPROX_MASK_COUNT * 9 * sizeof(int), NULL, &ret);
    cl_mem values_m = clCreateBuffer(program->context, CL_MEM_READ_ONLY, valueCount * 2 * vals_size, NULL, &ret);

    // write-only parameters
    cl_mem ac_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, valueCount * APPROX_MASK_COUNT * vals_size, NULL, &ret);
    cl_mem bc_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, valueCount * APPROX_MASK_COUNT * vals_size, NULL, &ret);
    cl_mem cc_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, valueCount * APPROX_MASK_COUNT * vals_size, NULL, &ret);

    // copy local data to GPU memory
    ret = clEnqueueWriteBuffer(program->commandQueue, mask_weights_m, CL_TRUE, 0, APPROX_MASK_COUNT * sizeof(int), mask_weights, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(program->commandQueue, mask_shift_base_m, CL_TRUE, 0, APPROX_MASK_COUNT * 9 * sizeof(int), mask_shift_base, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(program->commandQueue, values_m, CL_TRUE, 0, valueCount * 2 * vals_size, vals_cp, 0, NULL, NULL);

    cl_event event1;

    // create OpenCL kernel for derivatives calculation
    cl_kernel kernel = clCreateKernel(program->prog, "quadspline_calc_parameters", &ret);

    // set the arguments
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&mask_weights_m);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mask_shift_base_m);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&values_m);
    ret = clSetKernelArg(kernel, 3, sizeof(int), &valueCount);
    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&ac_m);
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&bc_m);
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&cc_m);

    // work!
    size_t local_item_size = 1;
    size_t global_item_size = APPROX_MASK_COUNT;

    ret = clEnqueueNDRangeKernel(program->commandQueue, kernel, 1, NULL, &global_item_size, &local_item_size, 0, NULL, &event1);

    clWaitForEvents(1, &event1);

    for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
    {
        aCoefs[mask].resize(valueCount);
        bCoefs[mask].resize(valueCount);
        cCoefs[mask].resize(valueCount);
    }

    // copy calculated parameters to inner vectors to make further value calculation possible
    // also distinguish by floating point number length (precision)
    if (clSupportsDouble())
    {
        double *tmp = new double[valueCount * APPROX_MASK_COUNT];
        ret = clEnqueueReadBuffer(program->commandQueue, ac_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                aCoefs[mask][j] = tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, bc_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                bCoefs[mask][j] = tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, cc_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                cCoefs[mask][j] = tmp[mask * valueCount + j];

        delete[] tmp;
    }
    else
    {
        float *tmp = new float[valueCount * APPROX_MASK_COUNT];
        ret = clEnqueueReadBuffer(program->commandQueue, ac_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                aCoefs[mask][j] = (double)tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, bc_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                bCoefs[mask][j] = (double)tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, cc_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                cCoefs[mask][j] = (double)tmp[mask * valueCount + j];

        delete[] tmp;
    }

    // flush and release local resources

    ret = clFlush(program->commandQueue);
    ret = clFinish(program->commandQueue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseMemObject(mask_weights_m);
    ret = clReleaseMemObject(mask_shift_base_m);
    ret = clReleaseMemObject(values_m);
    ret = clReleaseMemObject(ac_m);
    ret = clReleaseMemObject(bc_m);
    ret = clReleaseMemObject(cc_m);
}

HRESULT IfaceCalling ApproxQuadraticSpline::Approximate(TApproximationParams*)
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
    else if (appConcurrency == ConcurrencyType::ct_parallel_opencl)
    {
        CalculateParameters_OpenCL();
    }

    CLogical_Clock::Signal_Clock();

    return S_OK;
}

HRESULT IfaceCalling ApproxQuadraticSpline::GetBounds(TGlucoseLevelBounds *bounds)
{
    size_t remCount, maskedCount;

    remCount = (valueCount / 8);
    maskedCount = remCount * mask_weights[appCurrentTestMask];
    for (int i = 0; i < valueCount - remCount * 8; i++)
        maskedCount += (appCurrentTestMask >> (7 - i)) & 1;

    if (maskedCount < 1)
        return S_FALSE;

    size_t skipSides[2] = { 0, maskedCount };

    CFindMaskedBounds fb(values, appCurrentTestMask, skipSides);
    tbb::parallel_reduce(tbb::blocked_range<size_t>(0, maskedCount), fb);

    *bounds = fb.mBounds;

    return S_OK;
}

HRESULT IfaceCalling ApproxQuadraticSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder)
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
            if (index % mask_weights[mask] == 0 && index != 0)
                base += offsets[j + 1] - offsets[0];

            j = (j + 1) % mask_weights[mask];
        }
    }

    *filled = i;

    return S_OK;
}

void ApproxQuadraticSpline::CalculateCoefsFor(const uint32_t mask, size_t index, floattype aPrev, floattype bPrev, floattype xCur, floattype xNext, floattype yCur, floattype yNext)
{
    aCoefs[mask][index] = (((yCur - yNext) / (xCur - xNext)) - 2.0 * aPrev * xCur - bPrev) / (xNext - xCur);
    bCoefs[mask][index] = xCur * (2.0 * aPrev - 2.0 * aCoefs[mask][index]) + bPrev;
    cCoefs[mask][index] = yNext - aCoefs[mask][index] * xNext * xNext - bCoefs[mask][index] * xNext;
}

HRESULT ApproxQuadraticSpline::GetIndexFor(floattype time, size_t &index, size_t &origIndex, uint8_t mask)
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
