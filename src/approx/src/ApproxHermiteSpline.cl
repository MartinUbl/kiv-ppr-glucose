
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

__kernel void hermite_calc_parameters(__global const int *mask_weights, __global const int *mask_index_transform, __global const real_t *values, const int valueCount,
                                         __global real_t *xcoefs, __global real_t *ycoefs)
{
    int mask = get_global_id(0);
    int myWork = get_global_id(1);

    int pos0, pos1, pos2, pos3, base;

    pos3 = 8 * ((myWork + 3) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 3) % mask_weights[mask])];

    if (pos3 >= valueCount)
        return;

    pos2 = 8 * ((myWork + 2) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 2) % mask_weights[mask])];
    pos1 = 8 * ((myWork + 1) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 1) % mask_weights[mask])];
    pos0 = 8 * (myWork / mask_weights[mask]) + mask_index_transform[mask * 8 + (myWork % mask_weights[mask])];

    base = (mask * valueCount * 4) + (myWork + 1) * 4;

    xcoefs[base + 0] = values[pos1 * 2 + 0];
    xcoefs[base + 1] = -values[pos0 * 2 + 0] / 2.0 + values[pos2 * 2 + 0] / 2.0;
    xcoefs[base + 2] = values[pos0 * 2 + 0] - (5.0*values[pos1 * 2 + 0]) / 2.0 + 2.0*values[pos2 * 2 + 0] - values[pos3 * 2 + 0] / 2.0;
    xcoefs[base + 3] = -values[pos0 * 2 + 0] / 2.0 + (3.0f*values[pos1 * 2 + 0]) / 2.0 - (3.0*values[pos2 * 2 + 0]) / 2.0 + values[pos3 * 2 + 0] / 2.0;

    ycoefs[base + 0] = values[pos1 * 2 + 1];
    ycoefs[base + 1] = -values[pos0 * 2 + 1] / 2.0 + values[pos2 * 2 + 1] / 2.0;
    ycoefs[base + 2] = values[pos0 * 2 + 1] - (5.0*values[pos1 * 2 + 1]) / 2.0 + 2.0*values[pos2 * 2 + 1] - values[pos3 * 2 + 1] / 2.0f;
    ycoefs[base + 3] = -values[pos0 * 2 + 1] / 2.0 + (3.0*values[pos1 * 2 + 1]) / 2.0 - (3.0f*values[pos2 * 2 + 1]) / 2.0 + values[pos3 * 2 + 1] / 2.0;
}
