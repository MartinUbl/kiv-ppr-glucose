#pragma once

#include "CommonApprox.h"

#include <vector>

#define APPROX_MASK_COUNT 0x100

#pragma warning( push )
#pragma warning( disable : 4250 )   // C4250 - 'class1' : inherits 'class2::member' via dominance

class ApproxQuadraticSpline : public CCommonApprox
{
    public:
        ApproxQuadraticSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {};

        // Overriden methods from IApproximatedGlucoseLevels, see ApproxIface.h

        HRESULT IfaceCalling Approximate(TApproximationParams *params);
        HRESULT IfaceCalling GetLevels(floattype desiredtime, floattype stepping, size_t count,
                                       floattype *levels, size_t *filled, size_t derivationorder);

    protected:

        void CalculateParametersForMask(const uint32_t mask);

        void CalculateParameters_threads();

        void CalculateParameters_AMP();

        void CalculateParameters_TBB();

        void CalculateParameters_OpenCL();

        // calculates coefficients for one specific index in one specific mask
        inline void CalculateCoefsFor(const uint32_t mask, size_t index, floattype aPrev, floattype bPrev, floattype xCur, floattype xNext, floattype yCur, floattype yNext);

        // retrieves index of closest value in glucose levels
        HRESULT GetIndexFor(floattype time, size_t &index, size_t &origIndex, uint8_t mask = 0xFF);

        // quadratic spline coefficients
        std::vector<floattype> aCoefs[APPROX_MASK_COUNT], bCoefs[APPROX_MASK_COUNT], cCoefs[APPROX_MASK_COUNT];

        // cached value count
        size_t valueCount;

        // cached value pointer
        TGlucoseLevel* values;
};

#pragma warning( pop )
