
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

__kernel void quadspline_calc_parameters(__global const int *mask_weights, __global const int *mask_shift_base, __global const real_t *values, const int valueCount,
                                         __global real_t* aCoefs, __global real_t* bCoefs, __global real_t* cCoefs)
{
    int mask = get_global_id(0);
    int i, j = 1, base = 0;
    real_t aPrev, bPrev, xCur, xNext, yCur, yNext;

    aCoefs[mask * valueCount + 0] = 0;
    bCoefs[mask * valueCount + 0] = (values[2*mask_shift_base[mask * 9 + 1] + 1] - values[2*mask_shift_base[mask * 9 + 0] + 1]) / (values[2*mask_shift_base[mask * 9 + 1] + 0] - values[2*mask_shift_base[mask * 9 + 0] + 0]);
    cCoefs[mask * valueCount + 0] = values[2*mask_shift_base[mask * 9 + 0] + 1] - bCoefs[mask * valueCount + 0] * values[2*mask_shift_base[mask * 9 + 0] + 0];

    for (i = 1; i < valueCount - 1; i++, j++)
    {
        if (j == mask_weights[mask])
        {
            base += mask_shift_base[mask * 9 + j] - mask_shift_base[mask * 9 + 0];
            j = 0;
        }

        if (base + mask_shift_base[mask * 9 + j] > valueCount || base + mask_shift_base[mask * 9 + j + 1] > valueCount)
            break;

        aPrev = aCoefs[mask * valueCount + i - 1];
        bPrev = bCoefs[mask * valueCount + i - 1];

        xCur = values[2*(base + mask_shift_base[mask * 9 + j]) + 0];
        xNext = values[2*(base + mask_shift_base[mask * 9 + j + 1]) + 0];
        yCur = values[2*(base + mask_shift_base[mask * 9 + j]) + 1];
        yNext = values[2*(base + mask_shift_base[mask * 9 + j + 1]) + 1];

        aCoefs[mask * valueCount + i] = (((yCur - yNext) / (xCur - xNext)) - 2.0 * aPrev * xCur - bPrev) / (xNext - xCur);
        bCoefs[mask * valueCount + i] = xCur * (2.0 * aPrev - 2.0 * aCoefs[mask * valueCount + i]) + bPrev;
        cCoefs[mask * valueCount + i] = yNext - aCoefs[mask * valueCount + i] * xNext * xNext - bCoefs[mask * valueCount + i] * xNext;
    }
}
