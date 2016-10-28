
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

#define LINEAR_2D_INDEX(x,xdim,y) (x * xdim + y)
#define LINEAR_2D_INDEX_M(x,xdim,y) (x * xdim + y + 2)

__kernel void akima_calc_derivatives(__global const int *mask_weights, __global const int *mask_index_transform, __global real_t *m, __global real_t *values, int valueCount)
{
    int mask = get_global_id(0);
    int myWork = get_global_id(1);

    int pos, pos_n;

    pos = 8 * (myWork / mask_weights[mask]) + mask_index_transform[mask * 8 + (myWork % mask_weights[mask])];
    pos_n = 8 * ((myWork + 1) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 1) % mask_weights[mask])];

    if (pos >= valueCount || pos_n >= valueCount)
        return;

    m[LINEAR_2D_INDEX_M(mask, valueCount, myWork)] = (values[pos_n * 2 + 1] - values[pos * 2 + 1]) / (values[pos_n * 2 + 0] - values[pos * 2 + 0]);
}

__kernel void akima_fill_boundary_derivatives(__global float *m, __global int *maskedCounts, int valueCount)
{
    int mask = get_global_id(0);

    m[mask * valueCount + 1] = 2 * m[mask * valueCount + 2] - m[mask * valueCount + 3]; // "-1"
    m[mask * valueCount + 0] = 2 * m[mask * valueCount + 1] - m[mask * valueCount + 2]; // "-2"
}

__kernel void akima_calc_parameters(__global const int *mask_weights, __global const int *mask_index_transform, __global real_t *m, __global real_t *values, int valueCount,
                                    __global real_t *a0, __global real_t *a1, __global real_t *a2, __global real_t *a3)
{
    int mask = get_global_id(0);
    int myWork = get_global_id(1);

    int pos, pos_n;
    real_t tmpWeight1, tmpWeight2;

    pos = 8 * (myWork / mask_weights[mask]) + mask_index_transform[mask * 8 + (myWork % mask_weights[mask])];
    pos_n = 8 * ((myWork + 1) / mask_weights[mask]) + mask_index_transform[mask * 8 + ((myWork + 1) % mask_weights[mask])];

    if (pos >= valueCount || pos_n >= valueCount)
        return;

    real_t wsum1 = fabs(m[LINEAR_2D_INDEX_M(mask, valueCount, myWork + 1)] - m[LINEAR_2D_INDEX_M(mask, valueCount, myWork + 0)]) + fabs(m[LINEAR_2D_INDEX_M(mask, valueCount, myWork - 1)] - m[LINEAR_2D_INDEX_M(mask, valueCount, myWork - 2)]);
    if (wsum1 == 0)
    {
        a1[LINEAR_2D_INDEX(mask, valueCount, myWork)] = m[LINEAR_2D_INDEX_M(mask, valueCount, myWork)];
        a2[LINEAR_2D_INDEX(mask, valueCount, myWork)] = 0;
        a3[LINEAR_2D_INDEX(mask, valueCount, myWork)] = 0;
    }
    else
    {
        real_t deltax = values[pos_n * 2 + 0] - values[pos * 2 + 0];

        tmpWeight1 = fabs(m[LINEAR_2D_INDEX_M(mask, valueCount, myWork + 1)] - m[LINEAR_2D_INDEX_M(mask, valueCount, myWork + 0)]); // w_i+1
        tmpWeight2 = fabs(m[LINEAR_2D_INDEX_M(mask, valueCount, myWork - 1)] - m[LINEAR_2D_INDEX_M(mask, valueCount, myWork - 2)]); // w_i-1

        real_t wsum = tmpWeight1 + tmpWeight2;

        if (wsum == 0)
        {
            tmpWeight1 = 1.0;
            tmpWeight2 = 1.0;
            wsum = 2.0;
        }

        a0[LINEAR_2D_INDEX(mask, valueCount, myWork)] = values[pos * 2 + 1];
        a1[LINEAR_2D_INDEX(mask, valueCount, myWork)] = (tmpWeight1 * m[LINEAR_2D_INDEX_M(mask, valueCount, myWork - 1)] + tmpWeight2 * m[LINEAR_2D_INDEX_M(mask, valueCount, myWork)]) / wsum;

        tmpWeight1 = fabs(m[LINEAR_2D_INDEX_M(mask, valueCount, myWork + 1 + 1)] - m[LINEAR_2D_INDEX_M(mask, valueCount, myWork + 1)]);     // w_i+1
        tmpWeight2 = fabs(m[LINEAR_2D_INDEX_M(mask, valueCount, myWork - 1 + 1)] - m[LINEAR_2D_INDEX_M(mask, valueCount, myWork - 2 + 1)]); // w_i-1

        wsum = tmpWeight1 + tmpWeight2;

        if (wsum == 0)
        {
            tmpWeight1 = 1.0;
            tmpWeight2 = 1.0;
            wsum = 2.0;
        }

        real_t nextA1 = (tmpWeight1 * m[LINEAR_2D_INDEX_M(mask, valueCount, myWork - 1 + 1)] + tmpWeight2 * m[LINEAR_2D_INDEX_M(mask, valueCount, myWork + 1)]) / wsum;

        a2[LINEAR_2D_INDEX(mask, valueCount, myWork)] = (3 * m[LINEAR_2D_INDEX_M(mask, valueCount, myWork)] - 2 * a1[LINEAR_2D_INDEX(mask, valueCount, myWork)] - nextA1) / deltax;
        a3[LINEAR_2D_INDEX(mask, valueCount, myWork)] = (a1[LINEAR_2D_INDEX(mask, valueCount, myWork)] + nextA1 - 2 * m[LINEAR_2D_INDEX_M(mask, valueCount, myWork)]) / (deltax * deltax);
    }
}
