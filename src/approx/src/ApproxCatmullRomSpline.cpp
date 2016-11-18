#include "ApproxCatmullRomSpline.h"
#include "CommonApproxSupport.h"

#include <iostream>
#include <thread>
#include "../../cli/appconfig.h"
#include "precalculated.h"

#include <amp.h>
#include <tbb/tbb.h>

#include "../../cli/OpenCLLoader.h"

// safe, no name conflict
using namespace concurrency;

inline floattype getVectorSquareLength(const floattype x1, const floattype y1, const floattype x2, const floattype y2)
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

    const floattype tauCoef = tensionParameter;

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

        /*
        if (dt1 < FLT_EPSILON)
            dt1 = 1.0f;
        if (dt0 < FLT_EPSILON)
            dt0 = dt1;
        if (dt2 < FLT_EPSILON)
            dt2 = dt1;
        */

        // init xCoefs

        t1 = (values[pos1].datetime - values[pos0].datetime) / dt0 - (values[pos2].datetime - values[pos0].datetime) / (dt0 + dt1) + (values[pos2].datetime - values[pos1].datetime) / dt1;
        t2 = (values[pos2].datetime - values[pos1].datetime) / dt1 - (values[pos3].datetime - values[pos1].datetime) / (dt1 + dt2) + (values[pos3].datetime - values[pos2].datetime) / dt2;

        t1 *= dt1;
        t2 *= dt1;

        floattype* const xC = &xCoefs[mask][i + 1][0];

        xC[0] = values[pos1].datetime;
        xC[1] = t1;
        xC[2] = -3 * values[pos1].datetime + 3 * values[pos2].datetime - 2 * t1 - t2;
        xC[3] = 2 * values[pos1].datetime - 2 * values[pos2].datetime + t1 + t2;

        // init yCoefs

        t1 = (values[pos1].level - values[pos0].level) / dt0 - (values[pos2].level - values[pos0].level) / (dt0 + dt1) + (values[pos2].level - values[pos1].level) / dt1;
        t2 = (values[pos2].level - values[pos1].level) / dt1 - (values[pos3].level - values[pos1].level) / (dt1 + dt2) + (values[pos3].level - values[pos2].level) / dt2;

        t1 *= dt1;
        t2 *= dt1;

        floattype* const yC = &yCoefs[mask][i + 1][0];

        yC[0] = values[pos1].level;
        yC[1] = t1;
        yC[2] = -3 * values[pos1].level + 3 * values[pos2].level - 2 * t1 - t2;
        yC[3] = 2 * values[pos1].level - 2 * values[pos2].level + t1 + t2;
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

    size_t i;

    size_t maskedCount[APPROX_MASK_COUNT];
    size_t remCount;

    // serial calculation, this is too short for parallelization
    for (uint16_t mask = 1; mask < APPROX_MASK_COUNT; mask++)
    {
        remCount = (valueCount / 8);
        maskedCount[mask] = remCount * mask_weights[mask];
        for (i = 0; i < valueCount - remCount * 8; i++)
            maskedCount[mask] += (mask >> (7 - i)) & 1;

        xCoefs[mask].resize(maskedCount[mask] - 1);
        yCoefs[mask].resize(maskedCount[mask] - 1);
    }

    threadWork_t thrWork(maskedCount);

    for (i = 0; i < appWorkerCount; i++)
    {
        workers[i] = new std::thread([&]() {
            uint8_t mask;
            int32_t myWork;
            size_t pos0, pos1, pos2, pos3;
            floattype dt0, dt1, dt2, t1, t2;


            while (thrWork.getWork(mask, myWork))
            {
                pos3 = 8 * ((myWork + 3) / mask_weights[mask]) + mask_index_transform[mask][(myWork + 3) % mask_weights[mask]];

                if (pos3 >= valueCount)
                    continue;

                pos2 = 8 * ((myWork + 2) / mask_weights[mask]) + mask_index_transform[mask][(myWork + 2) % mask_weights[mask]];
                pos1 = 8 * ((myWork + 1) / mask_weights[mask]) + mask_index_transform[mask][(myWork + 1) % mask_weights[mask]];
                pos0 = 8 * (myWork / mask_weights[mask]) + mask_index_transform[mask][myWork % mask_weights[mask]];

                dt0 = pow(getVectorSquareLength(values[pos0].datetime, values[pos0].level, values[pos1].datetime, values[pos1].level), tensionParameter);
                dt1 = pow(getVectorSquareLength(values[pos1].datetime, values[pos1].level, values[pos2].datetime, values[pos2].level), tensionParameter);
                dt2 = pow(getVectorSquareLength(values[pos2].datetime, values[pos2].level, values[pos3].datetime, values[pos3].level), tensionParameter);

                /*
                if (dt1 < FLT_EPSILON)
                dt1 = 1.0f;
                if (dt0 < FLT_EPSILON)
                dt0 = dt1;
                if (dt2 < FLT_EPSILON)
                dt2 = dt1;
                */

                // init xCoefs

                t1 = (values[pos1].datetime - values[pos0].datetime) / dt0 - (values[pos2].datetime - values[pos0].datetime) / (dt0 + dt1) + (values[pos2].datetime - values[pos1].datetime) / dt1;
                t2 = (values[pos2].datetime - values[pos1].datetime) / dt1 - (values[pos3].datetime - values[pos1].datetime) / (dt1 + dt2) + (values[pos3].datetime - values[pos2].datetime) / dt2;

                t1 *= dt1;
                t2 *= dt1;

                floattype* const xC = &xCoefs[mask][myWork + 1][0];

                xC[0] = values[pos1].datetime;
                xC[1] = t1;
                xC[2] = -3 * values[pos1].datetime + 3 * values[pos2].datetime - 2 * t1 - t2;
                xC[3] = 2 * values[pos1].datetime - 2 * values[pos2].datetime + t1 + t2;

                // init yCoefs

                t1 = (values[pos1].level - values[pos0].level) / dt0 - (values[pos2].level - values[pos0].level) / (dt0 + dt1) + (values[pos2].level - values[pos1].level) / dt1;
                t2 = (values[pos2].level - values[pos1].level) / dt1 - (values[pos3].level - values[pos1].level) / (dt1 + dt2) + (values[pos3].level - values[pos2].level) / dt2;

                t1 *= dt1;
                t2 *= dt1;

                floattype* const yC = &yCoefs[mask][myWork + 1][0];

                yC[0] = values[pos1].level;
                yC[1] = t1;
                yC[2] = -3 * values[pos1].level + 3 * values[pos2].level - 2 * t1 - t2;
                yC[3] = 2 * values[pos1].level - 2 * values[pos2].level + t1 + t2;
            }
        });
    }

    for (i = 0; i < appWorkerCount; i++)
    {
        if (workers[i]->joinable())
            workers[i]->join();

        delete workers[i];
    }

    delete workers;
}

void ApproxCatmullRomSpline::CalculateParameters_AMP()
{
    //
}

void ApproxCatmullRomSpline::CalculateParameters_TBB()
{
    size_t maskedCount[APPROX_MASK_COUNT];

    // calculate masked counts
    tbb::parallel_for(size_t(1), size_t(APPROX_MASK_COUNT) - 1, size_t(1), [&](size_t mask) {
        size_t remCount = valueCount >> 3;// (valueCount / 8);
        maskedCount[mask] = remCount * mask_weights[mask];
        for (int i = 0; i < valueCount - remCount * 8; i++)
            maskedCount[mask] += (mask >> (7 - i)) & 1;
    });
    maskedCount[255] = valueCount;

    // resize value and derivation vectors, we know how much values do we need
    tbb::parallel_for(size_t(1), size_t(APPROX_MASK_COUNT), size_t(1), [&](size_t mask) {
        xCoefs[mask].resize(maskedCount[mask] - 1);
        yCoefs[mask].resize(maskedCount[mask] - 1);
    });

    // count derivations
    tbb::parallel_for(tbb::blocked_range2d<size_t, int>(1, APPROX_MASK_COUNT, 0, (int)valueCount), [&](const tbb::blocked_range2d<size_t, int> &idx) {
        size_t mask;
        int i;
        floattype dt0, dt1, dt2, t1, t2;
        size_t pos0, pos1, pos2, pos3;

        const int workbegin = idx.cols().begin();

        // for each mask
        for (mask = idx.rows().begin(); mask < idx.rows().end(); mask++)
        {
            const int limit = (int)min(idx.cols().end(), maskedCount[mask] - 1);

            // for each value
            for (i = workbegin; i < limit; i++)
            {
                pos3 = 8 * ((i + 3) / mask_weights[mask]) + mask_index_transform[mask][(i + 3) % mask_weights[mask]];

                if (pos3 >= valueCount)
                    continue;

                pos2 = 8 * ((i + 2) / mask_weights[mask]) + mask_index_transform[mask][(i + 2) % mask_weights[mask]];
                pos1 = 8 * ((i + 1) / mask_weights[mask]) + mask_index_transform[mask][(i + 1) % mask_weights[mask]];
                pos0 = 8 * (i / mask_weights[mask]) + mask_index_transform[mask][i % mask_weights[mask]];

                dt0 = pow(getVectorSquareLength(values[pos0].datetime, values[pos0].level, values[pos1].datetime, values[pos1].level), tensionParameter);
                dt1 = pow(getVectorSquareLength(values[pos1].datetime, values[pos1].level, values[pos2].datetime, values[pos2].level), tensionParameter);
                dt2 = pow(getVectorSquareLength(values[pos2].datetime, values[pos2].level, values[pos3].datetime, values[pos3].level), tensionParameter);

                /*
                if (dt1 < FLT_EPSILON)
                dt1 = 1.0f;
                if (dt0 < FLT_EPSILON)
                dt0 = dt1;
                if (dt2 < FLT_EPSILON)
                dt2 = dt1;
                */

                // init xCoefs

                t1 = (values[pos1].datetime - values[pos0].datetime) / dt0 - (values[pos2].datetime - values[pos0].datetime) / (dt0 + dt1) + (values[pos2].datetime - values[pos1].datetime) / dt1;
                t2 = (values[pos2].datetime - values[pos1].datetime) / dt1 - (values[pos3].datetime - values[pos1].datetime) / (dt1 + dt2) + (values[pos3].datetime - values[pos2].datetime) / dt2;

                t1 *= dt1;
                t2 *= dt1;

                floattype* const xC = &xCoefs[mask][i + 1][0];

                xC[0] = values[pos1].datetime;
                xC[1] = t1;
                xC[2] = -3 * values[pos1].datetime + 3 * values[pos2].datetime - 2 * t1 - t2;
                xC[3] = 2 * values[pos1].datetime - 2 * values[pos2].datetime + t1 + t2;

                // init yCoefs

                t1 = (values[pos1].level - values[pos0].level) / dt0 - (values[pos2].level - values[pos0].level) / (dt0 + dt1) + (values[pos2].level - values[pos1].level) / dt1;
                t2 = (values[pos2].level - values[pos1].level) / dt1 - (values[pos3].level - values[pos1].level) / (dt1 + dt2) + (values[pos3].level - values[pos2].level) / dt2;

                t1 *= dt1;
                t2 *= dt1;

                floattype* const yC = &yCoefs[mask][i + 1][0];

                yC[0] = values[pos1].level;
                yC[1] = t1;
                yC[2] = -3 * values[pos1].level + 3 * values[pos2].level - 2 * t1 - t2;
                yC[3] = 2 * values[pos1].level - 2 * values[pos2].level + t1 + t2;
            }
        }
    });
}

void ApproxCatmullRomSpline::CalculateParameters_OpenCL()
{
    cl_int ret;

    clProgramRecord* program = GetCLProgramRecord(apxmCatmullRomSpline);
    if (!program)
        return;

    int maskedCount[APPROX_MASK_COUNT];
    int remCount;

    maskedCount[0] = 0;
    // serial calculation, this is too short for parallelization
    for (uint16_t mask = 1; mask <= 255; mask++)
    {
        remCount = ((int)valueCount / 8);
        maskedCount[mask] = remCount * mask_weights[mask];
        for (int i = 0; i < valueCount - remCount * 8; i++)
            maskedCount[mask] += (mask >> (7 - i)) & 1;
    }

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
    cl_mem mask_index_transform_m = clCreateBuffer(program->context, CL_MEM_READ_ONLY, APPROX_MASK_COUNT * 8 * sizeof(int), NULL, &ret);
    cl_mem values_m = clCreateBuffer(program->context, CL_MEM_READ_ONLY, valueCount * 2 * vals_size, NULL, &ret);

    // write-only parameters
    cl_mem xcoefs_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, valueCount * APPROX_MASK_COUNT * vals_size * 4, NULL, &ret);
    cl_mem ycoefs_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, valueCount * APPROX_MASK_COUNT * vals_size * 4, NULL, &ret);

    // copy local data to GPU memory
    ret = clEnqueueWriteBuffer(program->commandQueue, mask_weights_m, CL_TRUE, 0, APPROX_MASK_COUNT * sizeof(int), mask_weights, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(program->commandQueue, mask_index_transform_m, CL_TRUE, 0, APPROX_MASK_COUNT * 8 * sizeof(int), mask_index_transform, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(program->commandQueue, values_m, CL_TRUE, 0, valueCount * 2 * vals_size, vals_cp, 0, NULL, NULL);

    cl_event event1;

    // create OpenCL kernel for derivatives calculation
    cl_kernel kernel = clCreateKernel(program->prog, "catmullrom_calc_parameters", &ret);

    // set the arguments
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&mask_weights_m);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mask_index_transform_m);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&values_m);
    ret = clSetKernelArg(kernel, 3, sizeof(int), &valueCount);
    ret = clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *)&xcoefs_m);
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&ycoefs_m);
    ret = clSetKernelArg(kernel, 6, sizeof(float), &tensionParameter);

    // work!
    size_t local_item_size[] = { 1, 16 };
    size_t global_item_size[] = { APPROX_MASK_COUNT, ShrRoundUp(16, valueCount) };

    ret = clEnqueueNDRangeKernel(program->commandQueue, kernel, 2, NULL, global_item_size, local_item_size, 0, NULL, &event1);

    clWaitForEvents(1, &event1);


    // resize destination vectors
    for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
    {
        xCoefs[mask].resize(maskedCount[mask]);
        yCoefs[mask].resize(maskedCount[mask]);
    }

    // copy calculated parameters to inner vectors to make further value calculation possible
    // also distinguish by floating point number length (precision)
    if (clSupportsDouble())
    {
        double *tmp = new double[valueCount * APPROX_MASK_COUNT * 4];
        ret = clEnqueueReadBuffer(program->commandQueue, xcoefs_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size * 4, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < maskedCount[mask]; j++)
                for (int k = 0; k < 4; k++)
                    xCoefs[mask][j][k] = tmp[mask * valueCount * 4 + j * 4 + k];

        ret = clEnqueueReadBuffer(program->commandQueue, ycoefs_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size * 4, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < maskedCount[mask]; j++)
                for (int k = 0; k < 4; k++)
                    yCoefs[mask][j][k] = tmp[mask * valueCount * 4 + j * 4 + k];

        free(tmp);
    }
    else
    {
        float *tmp = new float[valueCount * APPROX_MASK_COUNT * 4];
        ret = clEnqueueReadBuffer(program->commandQueue, xcoefs_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size * 4, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < maskedCount[mask]; j++)
                for (int k = 0; k < 4; k++)
                    xCoefs[mask][j][k] = (double)tmp[mask * valueCount * 4 + j * 4 + k];

        ret = clEnqueueReadBuffer(program->commandQueue, ycoefs_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size * 4, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < maskedCount[mask]; j++)
                for (int k = 0; k < 4; k++)
                    yCoefs[mask][j][k] = (double)tmp[mask * valueCount * 4 + j * 4 + k];

        free(tmp);
    }

    // flush and release local resources

    ret = clFlush(program->commandQueue);
    ret = clFinish(program->commandQueue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseMemObject(mask_weights_m);
    ret = clReleaseMemObject(mask_index_transform_m);
    ret = clReleaseMemObject(values_m);
    ret = clReleaseMemObject(xcoefs_m);
    ret = clReleaseMemObject(ycoefs_m);
}

HRESULT IfaceCalling ApproxCatmullRomSpline::Approximate(TApproximationParams *params)
{
    uint32_t mask;

    tensionParameter = params ? params->catmull.tensionParam : catmullRomParam_Centripetal;

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
        // for 250 values: mask-level: ~130ms, parameter-level: ~170ms on my developer machine (i5 4570)

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

        // we have to solve cubic equation for X coordinate in order to obtain real "time" value for Y coordinate calculation
        // this solution always lies between 0 and 1 (it's the "time" on curve, and since it's centripetal, it always has only one solution)
        calctime = Cubic_IdentitySolve(xCoefs[mask][index][3], xCoefs[mask][index][2], xCoefs[mask][index][1], xCoefs[mask][index][0] - curtime);

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
