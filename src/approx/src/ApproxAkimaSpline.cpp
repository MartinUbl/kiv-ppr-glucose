#include "ApproxAkimaSpline.h"

#include <iostream>
#include <thread>
#include "../../cli/appconfig.h"
#include "precalculated.h"

#include <amp.h>
#include <tbb/tbb.h>

#include "../../cli/OpenCLLoader.h"

// safe, no name conflict
using namespace concurrency;

void ApproxAkimaSpline::CalculateParametersForMask(const uint32_t mask)
{
    // shifting base for mask-based calculation
    size_t base;
    // offsets of shifted values (used for fast lookup) (one extra to avoid division by 8 in loop)
    const int *offsets = mask_shift_base[mask];

    size_t i, j;

    std::vector<floattype> derivations;

    size_t maskedCount;
    size_t remCount;

    floattype tmpWeight1, tmpWeight2;
    floattype* m;

    remCount = (valueCount / 8);
    maskedCount = remCount * mask_weights[mask];
    for (i = 0; i < valueCount - remCount * 8; i++)
        maskedCount += (mask >> (7 - i)) & 1;

    // resize value vectors, we know how much values do we need
    a0Coefs[mask].resize(maskedCount - 1);
    a1Coefs[mask].resize(maskedCount - 1);
    a2Coefs[mask].resize(maskedCount - 1);
    a3Coefs[mask].resize(maskedCount - 1);

    derivations.clear();
    derivations.resize(maskedCount + 4);

    m = &derivations[2];

    base = 0;
    j = 0;

    // precalculate derivatives
    for (i = 0; i < valueCount - 1; i++, j++)
    {
        if (j == mask_weights[mask])
        {
            base += offsets[j] - offsets[0];
            j = 0;
        }

        if (base + offsets[j] > maskedCount || base + offsets[j + 1] > maskedCount)
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
        if (j == mask_weights[mask])
        {
            base += offsets[j] - offsets[0];
            j = 0;
        }

        if (base + offsets[j] > maskedCount || base + offsets[j + 1] > maskedCount)
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

            tmpWeight1 = fabs(m[i + 1] - m[i]);     // w_i+1
            tmpWeight2 = fabs(m[i - 1] - m[i - 2]); // w_i-1

            floattype wsum = tmpWeight1 + tmpWeight2;

            if (wsum == 0)
            {
                tmpWeight1 = 1.0;
                tmpWeight2 = 1.0;
                wsum = 2.0;
            }

            a0Coefs[mask][i] = values[base + offsets[j]].level;
            a1Coefs[mask][i] = (tmpWeight1 * m[i - 1] + tmpWeight2 * m[i]) / wsum;

            tmpWeight1 = fabs(m[i + 1 + 1] - m[i + 1]);     // w_i+1
            tmpWeight2 = fabs(m[i - 1 + 1] - m[i - 2 + 1]); // w_i-1

            wsum = tmpWeight1 + tmpWeight2;

            if (wsum == 0)
            {
                tmpWeight1 = 1.0;
                tmpWeight2 = 1.0;
                wsum = 2.0;
            }

            const floattype nextA1 = (tmpWeight1 * m[i - 1 + 1] + tmpWeight2 * m[i + 1]) / wsum;

            a2Coefs[mask][i] = (3 * m[i] - 2 * a1Coefs[mask][i] - nextA1) / deltax;
            a3Coefs[mask][i] = (a1Coefs[mask][i] + nextA1 - 2 * m[i]) / (deltax * deltax);
        }
    }
}

void ApproxAkimaSpline::CalculateParameters_threads()
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

    std::vector<floattype> derivations[APPROX_MASK_COUNT];

    size_t maskedCount[APPROX_MASK_COUNT];
    size_t remCount;

    floattype* m[APPROX_MASK_COUNT];

    // serial calculation, this is too short for parallelization
    for (uint16_t mask = 1; mask < APPROX_MASK_COUNT; mask++)
    {
        remCount = (valueCount / 8);
        maskedCount[mask] = remCount * mask_weights[mask];
        for (i = 0; i < valueCount - remCount * 8; i++)
            maskedCount[mask] += (mask >> (7 - i)) & 1;

        // resize value vectors, we know how much values do we need
        a0Coefs[mask].resize(maskedCount[mask] - 1);
        a1Coefs[mask].resize(maskedCount[mask] - 1);
        a2Coefs[mask].resize(maskedCount[mask] - 1);
        a3Coefs[mask].resize(maskedCount[mask] - 1);

        derivations[mask].clear();
        derivations[mask].resize(maskedCount[mask] + 4);

        m[mask] = &derivations[mask][2];
    }

    threadWork_t thrWork(maskedCount);

    std::mutex barrierMtx;
    std::atomic<size_t> barrierCnt;
    std::condition_variable cvWorker, cvFarmer;

    int phase = 0;

    barrierCnt = appWorkerCount;

    for (i = 0; i < appWorkerCount; i++)
    {
        workers[i] = new std::thread([&](ApproxAkimaSpline* obj) {
            uint8_t mask;
            int32_t myWork;
            int pos, pos_n;
            while (thrWork.getWork(mask, myWork))
            {
                pos = 8 * (myWork / mask_weights[mask]) + mask_index_transform[mask][myWork % mask_weights[mask]];
                pos_n = 8 * ((myWork + 1) / mask_weights[mask]) + mask_index_transform[mask][(myWork + 1) % mask_weights[mask]];

                if (pos >= valueCount || pos_n >= valueCount)
                    continue;

                m[mask][myWork] = (obj->values[pos_n].level - obj->values[pos].level) /
                                  (obj->values[pos_n].datetime - obj->values[pos].datetime);
            }

            // block on condition variable
            {
                std::unique_lock<std::mutex> lck(barrierMtx);
                if (--barrierCnt == 0)
                    cvFarmer.notify_all();

                while (phase == 0)
                    cvWorker.wait(lck);
            }

            floattype tmpWeight1, tmpWeight2;
            while (thrWork.getWork(mask, myWork))
            {
                pos = 8 * (myWork / mask_weights[mask]) + mask_index_transform[mask][myWork % mask_weights[mask]];
                pos_n = 8 * ((myWork + 1) / mask_weights[mask]) + mask_index_transform[mask][(myWork + 1) % mask_weights[mask]];

                if (pos >= valueCount || pos_n >= valueCount)
                    continue;

                floattype* mmask = m[mask];

                floattype wsum1 = fabs(mmask[myWork + 1] - mmask[myWork]) + fabs(mmask[myWork - 1] - mmask[myWork - 2]);
                if (wsum1 == 0)
                {
                    obj->a1Coefs[mask][myWork] = mmask[myWork];
                    obj->a2Coefs[mask][myWork] = 0;
                    obj->a3Coefs[mask][myWork] = 0;
                }
                else
                {
                    floattype deltax = obj->values[pos_n].datetime - obj->values[pos].datetime;

                    tmpWeight1 = fabs(mmask[myWork + 1] - mmask[myWork]);     // w_i+1
                    tmpWeight2 = fabs(mmask[myWork - 1] - mmask[myWork - 2]); // w_i-1

                    floattype wsum = tmpWeight1 + tmpWeight2;

                    if (wsum == 0)
                    {
                        tmpWeight1 = 1.0;
                        tmpWeight2 = 1.0;
                        wsum = 2.0;
                    }

                    obj->a0Coefs[mask][myWork] = obj->values[pos].level;
                    obj->a1Coefs[mask][myWork] = (tmpWeight1 * mmask[myWork - 1] + tmpWeight2 * mmask[myWork]) / wsum;

                    tmpWeight1 = fabs(mmask[myWork + 1 + 1] - mmask[myWork + 1]);     // w_i+1
                    tmpWeight2 = fabs(mmask[myWork - 1 + 1] - mmask[myWork - 2 + 1]); // w_i-1

                    wsum = tmpWeight1 + tmpWeight2;

                    if (wsum == 0)
                    {
                        tmpWeight1 = 1.0;
                        tmpWeight2 = 1.0;
                        wsum = 2.0;
                    }

                    const floattype nextA1 = (tmpWeight1 * mmask[myWork - 1 + 1] + tmpWeight2 * mmask[myWork + 1]) / wsum;

                    obj->a2Coefs[mask][myWork] = (3 * mmask[myWork] - 2 * obj->a1Coefs[mask][myWork] - nextA1) / deltax;
                    obj->a3Coefs[mask][myWork] = (obj->a1Coefs[mask][myWork] + nextA1 - 2 * mmask[myWork]) / (deltax * deltax);
                }
            }
        }, this);
    }

    // barrier
    {
        std::unique_lock<std::mutex> lck(barrierMtx);
        while (barrierCnt > 0)
            cvFarmer.wait(lck);

        // calculate boundary derivatives
        for (uint16_t mask = 1; mask < APPROX_MASK_COUNT; mask++)
        {
            floattype* mmask = m[mask];
            mmask[-1] = 2 * mmask[0] - mmask[1];
            mmask[-2] = 2 * mmask[-1] - mmask[0];
            mmask[maskedCount[mask] - 1] = 2 * mmask[maskedCount[mask] - 2] - mmask[maskedCount[mask] - 3];
            mmask[maskedCount[mask]] = 2 * mmask[maskedCount[mask] - 1] - mmask[maskedCount[mask] - 2];
        }

        thrWork.resetWork();
        phase = 1;

        cvWorker.notify_all();
    }

    for (i = 0; i < appWorkerCount; i++)
    {
        if (workers[i]->joinable())
            workers[i]->join();

        delete workers[i];
    }

    delete[] workers;
}

void ApproxAkimaSpline::CalculateParameters_AMP()
{
    //
}

void ApproxAkimaSpline::CalculateParameters_TBB()
{
    int i;

    std::vector<floattype> derivations[APPROX_MASK_COUNT];

    size_t maskedCount[APPROX_MASK_COUNT];
    size_t remCount;

    floattype* m[APPROX_MASK_COUNT];

    // serial calculation, this is too short for parallelization
    for (uint16_t mask = 1; mask <= 255; mask++)
    {
        remCount = (valueCount / 8);
        maskedCount[mask] = remCount * mask_weights[mask];
        for (i = 0; i < valueCount - remCount * 8; i++)
            maskedCount[mask] += (mask >> (7 - i)) & 1;

        // resize value vectors, we know how much values do we need
        a0Coefs[mask].resize(maskedCount[mask] - 1);
        a1Coefs[mask].resize(maskedCount[mask] - 1);
        a2Coefs[mask].resize(maskedCount[mask] - 1);
        a3Coefs[mask].resize(maskedCount[mask] - 1);

        derivations[mask].clear();
        derivations[mask].resize(maskedCount[mask] + 4);

        m[mask] = &derivations[mask][2];
    }

    // count derivations
    tbb::parallel_for(tbb::blocked_range2d<size_t, int>(1, APPROX_MASK_COUNT, 0, (int)valueCount), [&](const tbb::blocked_range2d<size_t, int> &idx) {
        size_t mask;
        int i, pos, pos_n;

        // for each mask
        for (mask = idx.rows().begin(); mask < idx.rows().end(); mask++)
        {
            // for each value
            for (i = idx.cols().begin(); i < idx.cols().end(); i++)
            {
                pos = 8 * (i / mask_weights[mask]) + mask_index_transform[mask][i % mask_weights[mask]];
                pos_n = 8 * ((i + 1) / mask_weights[mask]) + mask_index_transform[mask][(i + 1) % mask_weights[mask]];

                if (pos >= valueCount || pos_n >= valueCount)
                    continue;

                m[mask][i] = (values[pos_n].level - values[pos].level) /
                             (values[pos_n].datetime - values[pos].datetime);
            }
        }
    });

    // calculate boundary derivatives
    for (uint16_t mask = 1; mask < APPROX_MASK_COUNT; mask++)
    {
        floattype* mmask = m[mask];
        mmask[-1] = 2 * mmask[0] - mmask[1];
        mmask[-2] = 2 * mmask[-1] - mmask[0];
        mmask[maskedCount[mask] - 1] = 2 * mmask[maskedCount[mask] - 2] - mmask[maskedCount[mask] - 3];
        mmask[maskedCount[mask]] = 2 * mmask[maskedCount[mask] - 1] - mmask[maskedCount[mask] - 2];
    }

    // count parameters
    tbb::parallel_for(tbb::blocked_range2d<size_t, int>(1, APPROX_MASK_COUNT, 0, (int)valueCount), [&](const tbb::blocked_range2d<size_t, int> &idx) {
        size_t mask;
        volatile int i, pos, pos_n;
        floattype tmpWeight1, tmpWeight2;

        // for each mask
        for (mask = idx.rows().begin(); mask < idx.rows().end(); mask++)
        {
            // for each value
            for (i = idx.cols().begin(); i < idx.cols().end(); i++)
            {
                pos = 8 * (i / mask_weights[mask]) + mask_index_transform[mask][i % mask_weights[mask]];
                pos_n = 8 * ((i + 1) / mask_weights[mask]) + mask_index_transform[mask][(i + 1) % mask_weights[mask]];

                if (pos >= valueCount || pos_n >= valueCount)
                    continue;

                floattype* mmask = m[mask];

                floattype wsum1 = fabs(mmask[i + 1] - mmask[i]) + fabs(mmask[i - 1] - mmask[i - 2]);
                if (wsum1 == 0)
                {
                    a1Coefs[mask][i] = mmask[i];
                    a2Coefs[mask][i] = 0;
                    a3Coefs[mask][i] = 0;
                }
                else
                {
                    floattype deltax = values[pos_n].datetime - values[pos].datetime;

                    tmpWeight1 = fabs(mmask[i + 1] - mmask[i]);     // w_i+1
                    tmpWeight2 = fabs(mmask[i - 1] - mmask[i - 2]); // w_i-1

                    floattype wsum = tmpWeight1 + tmpWeight2;

                    if (wsum == 0)
                    {
                        tmpWeight1 = 1.0;
                        tmpWeight2 = 1.0;
                        wsum = 2.0;
                    }

                    a0Coefs[mask][i] = values[pos].level;
                    a1Coefs[mask][i] = (tmpWeight1 * mmask[i - 1] + tmpWeight2 * mmask[i]) / wsum;

                    tmpWeight1 = fabs(mmask[i + 1 + 1] - mmask[i + 1]);     // w_i+1
                    tmpWeight2 = fabs(mmask[i - 1 + 1] - mmask[i - 2 + 1]); // w_i-1

                    wsum = tmpWeight1 + tmpWeight2;

                    if (wsum == 0)
                    {
                        tmpWeight1 = 1.0;
                        tmpWeight2 = 1.0;
                        wsum = 2.0;
                    }

                    const floattype nextA1 = (tmpWeight1 * mmask[i - 1 + 1] + tmpWeight2 * mmask[i + 1]) / wsum;

                    a2Coefs[mask][i] = (3 * mmask[i] - 2 * a1Coefs[mask][i] - nextA1) / deltax;
                    a3Coefs[mask][i] = (a1Coefs[mask][i] + nextA1 - 2 * mmask[i]) / (deltax * deltax);
                }
            }
        }
    });
}

void ApproxAkimaSpline::CalculateParameters_OpenCL()
{
    cl_int ret;

    clProgramRecord* program = GetCLProgramRecord(apxmAkimaSpline);
    if (!program)
        return;

    int maskedCount[APPROX_MASK_COUNT];
    int remCount;

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
    cl_mem masked_counts_m = clCreateBuffer(program->context, CL_MEM_READ_ONLY, APPROX_MASK_COUNT * sizeof(int), NULL, &ret);

    // write-only parameters
    cl_mem ders_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, (valueCount + 4) * APPROX_MASK_COUNT * vals_size, NULL, &ret);
    cl_mem a0_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, valueCount * APPROX_MASK_COUNT * vals_size, NULL, &ret);
    cl_mem a1_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, valueCount * APPROX_MASK_COUNT * vals_size, NULL, &ret);
    cl_mem a2_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, valueCount * APPROX_MASK_COUNT * vals_size, NULL, &ret);
    cl_mem a3_m = clCreateBuffer(program->context, CL_MEM_WRITE_ONLY, valueCount * APPROX_MASK_COUNT * vals_size, NULL, &ret);

    // copy local data to GPU memory
    ret = clEnqueueWriteBuffer(program->commandQueue, mask_weights_m, CL_TRUE, 0, APPROX_MASK_COUNT * sizeof(int), mask_weights, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(program->commandQueue, mask_index_transform_m, CL_TRUE, 0, APPROX_MASK_COUNT * 8 * sizeof(int), mask_index_transform, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(program->commandQueue, values_m, CL_TRUE, 0, valueCount * 2 * vals_size, vals_cp, 0, NULL, NULL);
    ret = clEnqueueWriteBuffer(program->commandQueue, masked_counts_m, CL_TRUE, 0, APPROX_MASK_COUNT, maskedCount, 0, NULL, NULL);

    cl_event event1, event2, event3;

    ///////////////////////////// derivatives calculation /////////////////////////////

    // create OpenCL kernel for derivatives calculation
    cl_kernel kernel = clCreateKernel(program->prog, "akima_calc_derivatives", &ret);

    // set the arguments
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&mask_weights_m);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mask_index_transform_m);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&ders_m);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&values_m);
    ret = clSetKernelArg(kernel, 4, sizeof(int), &valueCount);

    // work!
    size_t local_item_size[] = { 1, 16 };
    size_t global_item_size[] = { APPROX_MASK_COUNT, ShrRoundUp(16, valueCount) };

    ret = clEnqueueNDRangeKernel(program->commandQueue, kernel, 2, NULL, global_item_size, local_item_size, 0, NULL, &event1);

    clWaitForEvents(1, &event1);

    ////////////////////// derivatives boundary values calculation ////////////////////

    ret = clReleaseKernel(kernel);

    kernel = clCreateKernel(program->prog, "akima_fill_boundary_derivatives", &ret);
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&ders_m);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&masked_counts_m);
    ret = clSetKernelArg(kernel, 2, sizeof(int), &valueCount);

    size_t local_item_size2 = 16;
    size_t global_item_size2 = APPROX_MASK_COUNT;
    ret = clEnqueueNDRangeKernel(program->commandQueue, kernel, 1, NULL, &global_item_size2, &local_item_size2, 0, NULL, &event2);

    clWaitForEvents(1, &event2);

    ///////////////////////////// parameters calculation //////////////////////////////

    ret = clReleaseKernel(kernel);

    kernel = clCreateKernel(program->prog, "akima_calc_parameters", &ret);
    ret = clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *)&mask_weights_m);
    ret = clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *)&mask_index_transform_m);
    ret = clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *)&ders_m);
    ret = clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *)&values_m);
    ret = clSetKernelArg(kernel, 4, sizeof(int), &valueCount);
    ret = clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *)&a0_m);
    ret = clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *)&a1_m);
    ret = clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *)&a2_m);
    ret = clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *)&a3_m);

    size_t local_item_size3[] = { 1, 16 };
    size_t global_item_size3[] = { APPROX_MASK_COUNT, ShrRoundUp(16, valueCount) };
    ret = clEnqueueNDRangeKernel(program->commandQueue, kernel, 2, NULL, global_item_size3, local_item_size3, 0, NULL, &event3);

    clWaitForEvents(1, &event3);


    // resize destination vectors
    for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
    {
        a0Coefs[mask].resize(valueCount);
        a1Coefs[mask].resize(valueCount);
        a2Coefs[mask].resize(valueCount);
        a3Coefs[mask].resize(valueCount);
    }

    // copy calculated parameters to inner vectors to make further value calculation possible
    // also distinguish by floating point number length (precision)
    if (clSupportsDouble())
    {
        double *tmp = new double[valueCount * APPROX_MASK_COUNT];
        ret = clEnqueueReadBuffer(program->commandQueue, a0_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                a0Coefs[mask][j] = tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, a1_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                a1Coefs[mask][j] = tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, a2_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                a2Coefs[mask][j] = tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, a3_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                a3Coefs[mask][j] = tmp[mask * valueCount + j];

        free(tmp);
    }
    else
    {
        float *tmp = new float[valueCount * APPROX_MASK_COUNT];
        ret = clEnqueueReadBuffer(program->commandQueue, a0_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                a0Coefs[mask][j] = (double)tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, a1_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                a1Coefs[mask][j] = (double)tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, a2_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                a2Coefs[mask][j] = (double)tmp[mask * valueCount + j];

        ret = clEnqueueReadBuffer(program->commandQueue, a3_m, CL_TRUE, 0, valueCount * APPROX_MASK_COUNT * vals_size, tmp, 0, NULL, NULL);

        for (int mask = 0; mask < APPROX_MASK_COUNT; mask++)
            for (int j = 0; j < valueCount; j++)
                a3Coefs[mask][j] = (double)tmp[mask * valueCount + j];

        free(tmp);
    }

    // flush and release local resources

    ret = clFlush(program->commandQueue);
    ret = clFinish(program->commandQueue);
    ret = clReleaseKernel(kernel);
    ret = clReleaseMemObject(mask_weights_m);
    ret = clReleaseMemObject(mask_index_transform_m);
    ret = clReleaseMemObject(values_m);
    ret = clReleaseMemObject(masked_counts_m);
    ret = clReleaseMemObject(ders_m);
}

HRESULT IfaceCalling ApproxAkimaSpline::Approximate(TApproximationParams *params)
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

    return S_OK;
}

HRESULT IfaceCalling ApproxAkimaSpline::GetLevels(floattype desiredtime, floattype stepping, size_t count, floattype *levels, size_t *filled, size_t derivationorder)
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
    floattype timediff;

    // fill the array
    for (; i < count; i++)
    {
        // we want more than we have
        if (index >= maskedCount - 1)
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
            if (index % mask_weights[mask] == 0 && index != 0)
                base += offsets[j + 1] - offsets[0];

            j = (j + 1) % mask_weights[mask];
        }
    }

    *filled = i;

    return S_OK;
}

HRESULT ApproxAkimaSpline::GetIndexFor(floattype time, size_t &index, size_t &origIndex, uint8_t mask)
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
