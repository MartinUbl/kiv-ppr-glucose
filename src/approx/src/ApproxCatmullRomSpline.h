#pragma once

#include "CommonApprox.h"

#include <array>
#include <vector>
#include <thread>

const floattype catmullRomParam_Centripetal = 0.25;

#define APPROX_MASK_COUNT 0x100

#pragma warning( push )
#pragma warning( disable : 4250 )   // C4250 - 'class1' : inherits 'class2::member' via dominance

class ApproxCatmullRomSpline : public CCommonApprox
{
    public:
        ApproxCatmullRomSpline(IGlucoseLevels *levels) : CCommonApprox(levels) {};

        // Overriden methods from IApproximatedGlucoseLevels, see ApproxIface.h

        HRESULT IfaceCalling Approximate(TApproximationParams *params);
        HRESULT IfaceCalling GetLevels(floattype desiredtime, floattype stepping, size_t count,
                                       floattype *levels, size_t *filled, size_t derivationorder);
        HRESULT IfaceCalling GetBounds(TGlucoseLevelBounds *bounds);

    protected:

        // calculate parameters for single mask (serial version)
        void CalculateParametersForMask(const uint32_t mask);

        // calculate parameters using standard threads
        void CalculateParameters_threads();

        // calculate paàameters using C++AMP
        void CalculateParameters_AMP();

        // calculate parameters using Intel TBB
        void CalculateParameters_TBB();

        // calculate parameters using OpenCL
        void CalculateParameters_OpenCL();

        // retrieves index of closest value in glucose levels
        HRESULT GetIndexFor(floattype time, size_t &index, size_t &origIndex, uint8_t mask = 0xFF);

        // cattmul-rom spline coefficients
        std::vector<std::array<floattype, 4>> xCoefs[APPROX_MASK_COUNT], yCoefs[APPROX_MASK_COUNT];

        // cached value count
        size_t valueCount;

        // cached value pointer
        TGlucoseLevel* values;

        // Catmull-Rom polynomial formula tension parameter
        floattype tensionParameter;
};

#pragma warning( pop )
