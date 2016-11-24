
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
                                         __global real_t *xcoefs, __global real_t *ycoefs, const real_t tensionParam)
{
    int mask = get_global_id(0);
    int myWork = get_global_id(1);

    int pos0, pos1, pos2, pos3, base;
    real_t t1, t2, p1, p2;

    pos3 = 8 * ((myWork + 3) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 3) % mask_weights[mask])];

    if (pos3 >= valueCount)
        return;

    pos2 = 8 * ((myWork + 2) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 2) % mask_weights[mask])];
    pos1 = 8 * ((myWork + 1) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 1) % mask_weights[mask])];
    pos0 = 8 * (myWork / mask_weights[mask]) + mask_index_transform[mask * 8 + (myWork % mask_weights[mask])];

    base = (mask * valueCount * 4) + (myWork + 1) * 4;

    p1 = values[pos1 * 2 + 0];
    p2 = values[pos2 * 2 + 0];
    // tangent estimation using neighbors
    t1 = (1.0 - tensionParam) * (values[pos2 * 2 + 0] - values[pos0 * 2 + 0]);
    t2 = (1.0 - tensionParam) * (values[pos3 * 2 + 0] - values[pos1 * 2 + 0]);

    xcoefs[base + 3] = 2 * p1 - 2 * p2 + t1 + t2;
    xcoefs[base + 2] = -3 * p1 + 3 * p2 - 2 * t1 - t2;
    xcoefs[base + 1] = t1;
    xcoefs[base + 0] = p1;

    p1 = values[pos1 * 2 + 1];
    p2 = values[pos2 * 2 + 1];
    // tangent estimation using neighbors
    t1 = (1.0 - tensionParam) * (values[pos2 * 2 + 1] - values[pos0 * 2 + 1]);
    t2 = (1.0 - tensionParam) * (values[pos3 * 2 + 1] - values[pos1 * 2 + 1]);

    ycoefs[base + 3] = 2 * p1 - 2 * p2 + t1 + t2;
    ycoefs[base + 2] = -3 * p1 + 3 * p2 - 2 * t1 - t2;
    ycoefs[base + 1] = t1;
    ycoefs[base + 0] = p1;
}
