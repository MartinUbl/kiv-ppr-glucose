
#if defined(cl_khr_fp64)
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#define DOUBLE_SUPPORT_AVAILABLE
#endif

#if defined(DOUBLE_SUPPORT_AVAILABLE)
typedef double real_t;
#else
typedef float real_t;
#endif

real_t getVectorSquareLength(real_t x1, real_t y1, real_t x2, real_t y2)
{
    real_t dx = x2 - x1;
    real_t dy = y2 - y1;
    return dx*dx + dy*dy;
}

__kernel void catmullrom_calc_parameters(__global const int *mask_weights, __global const int *mask_index_transform, __global const real_t *values, const int valueCount,
                                         __global real_t *xcoefs, __global real_t *ycoefs, const float tensionParameter)
{
    int mask = get_global_id(0);
    int myWork = get_global_id(1);

    int pos0, pos1, pos2, pos3, base;
    real_t dt0, dt1, dt2, t1, t2;

    pos3 = 8 * ((myWork + 3) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 3) % mask_weights[mask])];

    if (pos3 >= valueCount)
        return;

    pos2 = 8 * ((myWork + 2) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 2) % mask_weights[mask])];
    pos1 = 8 * ((myWork + 1) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 1) % mask_weights[mask])];
    pos0 = 8 * (myWork / mask_weights[mask]) + mask_index_transform[mask * 8 + (myWork % mask_weights[mask])];

    dt0 = pow((real_t)getVectorSquareLength(values[pos0 * 2 + 0], values[pos0 * 2 + 1], values[pos1 * 2 + 0], values[pos1 * 2 + 1]), (real_t)tensionParameter);
    dt1 = pow((real_t)getVectorSquareLength(values[pos1 * 2 + 0], values[pos1 * 2 + 1], values[pos2 * 2 + 0], values[pos2 * 2 + 1]), (real_t)tensionParameter);
    dt2 = pow((real_t)getVectorSquareLength(values[pos2 * 2 + 0], values[pos2 * 2 + 1], values[pos3 * 2 + 0], values[pos3 * 2 + 1]), (real_t)tensionParameter);

    // init xCoefs

    t1 = (values[pos1 * 2 + 0] - values[pos0 * 2 + 0]) / dt0 - (values[pos2 * 2 + 0] - values[pos0 * 2 + 0]) / (dt0 + dt1) + (values[pos2 * 2 + 0] - values[pos1 * 2 + 0]) / dt1;
    t2 = (values[pos2 * 2 + 0] - values[pos1 * 2 + 0]) / dt1 - (values[pos3 * 2 + 0] - values[pos1 * 2 + 0]) / (dt1 + dt2) + (values[pos3 * 2 + 0] - values[pos2 * 2 + 0]) / dt2;

    t1 *= dt1;
    t2 *= dt1;

    base = (mask * valueCount * 4) + (myWork + 1) * 4;

    xcoefs[base + 0] = values[pos1 * 2 + 0];
    xcoefs[base + 1] = t1;
    xcoefs[base + 2] = -3 * values[pos1 * 2 + 0] + 3 * values[pos2 * 2 + 0] - 2 * t1 - t2;
    xcoefs[base + 3] = 2 * values[pos1 * 2 + 0] - 2 * values[pos2 * 2 + 0] + t1 + t2;

    // init yCoefs

    t1 = (values[pos1 * 2 + 1] - values[pos0 * 2 + 1]) / dt0 - (values[pos2 * 2 + 1] - values[pos0 * 2 + 1]) / (dt0 + dt1) + (values[pos2 * 2 + 1] - values[pos1 * 2 + 1]) / dt1;
    t2 = (values[pos2 * 2 + 1] - values[pos1 * 2 + 1]) / dt1 - (values[pos3 * 2 + 1] - values[pos1 * 2 + 1]) / (dt1 + dt2) + (values[pos3 * 2 + 1] - values[pos2 * 2 + 1]) / dt2;

    t1 *= dt1;
    t2 *= dt1;

    ycoefs[base + 0] = values[pos1 * 2 + 1];
    ycoefs[base + 1] = t1;
    ycoefs[base + 2] = -3 * values[pos1 * 2 + 1] + 3 * values[pos2 * 2 + 1] - 2 * t1 - t2;
    ycoefs[base + 3] = 2 * values[pos1 * 2 + 1] - 2 * values[pos2 * 2 + 1] + t1 + t2;
}
