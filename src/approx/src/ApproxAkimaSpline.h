#pragma once

#include "CommonApprox.h"

#include <vector>
#include <thread>

#define APPROX_MASK_COUNT 0x100

#pragma warning( push )
#pragma warning( disable : 4250 )   // C4250 - 'class1' : inherits 'class2::member' via dominance

class ApproxAkimaSpline : public CCommonApprox
{
    public:
        ApproxAkimaSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {};

        // Overriden methods from IApproximatedGlucoseLevels, see ApproxIface.h

        HRESULT IfaceCalling Approximate(TApproximationParams *params);
        HRESULT IfaceCalling GetLevels(floattype desiredtime, floattype stepping, size_t count,
                                       floattype *levels, size_t *filled, size_t derivationorder);

    protected:

        void CalculateParametersForMask(const uint32_t mask);

        void CalculateParameters_threads();

        void CalculateParameters_AMP();

        void CalculateParameters_TBB();

        // retrieves index of closest value in glucose levels
        HRESULT GetIndexFor(floattype time, size_t &index);

        // quadratic spline coefficients
        std::vector<floattype> a0Coefs[APPROX_MASK_COUNT], a1Coefs[APPROX_MASK_COUNT], a2Coefs[APPROX_MASK_COUNT], a3Coefs[APPROX_MASK_COUNT];

        // cached value count
        size_t valueCount;

        // cached value pointer
        TGlucoseLevel* values;
};

#pragma warning( pop )
