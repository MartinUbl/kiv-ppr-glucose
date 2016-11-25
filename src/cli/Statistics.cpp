#include "Statistics.h"
#include "appconfig.h"
#include "../approx/src/precalculated.h"
#include <iostream>
#include <algorithm>

enum ErrorClassType
{
    ERR_ALL = 0,
    ERR_ZERO = 1,
    ERR_ONE = 2
};

void CalculateStats(std::vector<floattype> &errors, floattype &avg, floattype &minval, floattype &fquart, floattype &med, floattype &tquart, floattype &maxval, floattype &deviation)
{
    const size_t &levCount = errors.size();

    // may happen in mask 255 in zero masked values
    if (levCount == 0)
    {
        avg = 0;
        minval = 0;
        maxval = 0;
        fquart = 0;
        med = 0;
        tquart = 0;
        deviation = 0;
        return;
    }

    std::sort(errors.begin(), errors.end());

    minval = errors[0];
    maxval = errors[levCount - 1];

    if (levCount % 2 == 0)
        med = (errors[levCount / 2] + errors[1 + (levCount / 2)]) / 2;
    else
        med = errors[levCount / 2];

    if (levCount % 4 == 0)
    {
        fquart = (errors[levCount / 4] + errors[1 + (levCount / 4)]) / 2;
        tquart = (errors[3 * levCount / 4] + errors[1 + (3 * levCount / 4)]) / 2;
    }
    else
    {
        fquart = errors[levCount / 4];
        tquart = errors[3 * levCount / 4];
    }

    avg = 0;
    for (int i = 0; i < levCount; i++)
        avg += errors[i];
    avg /= (floattype)levCount;

    deviation = 0;
    for (int i = 0; i < levCount; i++)
        deviation += pow(errors[i] - avg, 2);
    deviation = sqrt(deviation / (floattype)levCount);
}

void CalculateAndPrintStatsForMask(size_t cur, std::vector<CGlucoseLevels*> &levelsVector, std::vector<CCommonApprox*> &approxVector)
{
    int i;
    TGlucoseLevel* levs;
    size_t levCount, remCount, maskedCount, maskedZeroCount;
    size_t filled;
    floattype target, curtime;
    std::vector<floattype> errors[3], relErrors[3], errorsDer[3];
    std::vector<floattype> refDerivations;

    levelsVector[cur]->GetLevels(&levs);
    levelsVector[cur]->GetLevelsCount(&levCount);

    remCount = levCount >> 3; //(valueCount / 8);
    maskedCount = remCount * mask_weights[appCurrentTestMask];
    for (i = 0; i < levCount - remCount * 8; i++)
        maskedCount += (appCurrentTestMask >> (7 - i)) & 1;

    maskedZeroCount = levCount - maskedCount;

    refDerivations.resize(levCount);
    for (i = 0; i < levCount; i++)
    {
        curtime = levs[i].datetime;
        approxVector[cur]->GetLevels(curtime, 0, 1, &refDerivations[i], &filled, 1);
    }

    errors[ERR_ALL].resize(levCount);
    relErrors[ERR_ALL].resize(levCount);
    errorsDer[ERR_ALL].resize(levCount);
    errors[ERR_ZERO].resize(maskedZeroCount);
    relErrors[ERR_ZERO].resize(maskedZeroCount);
    errorsDer[ERR_ZERO].resize(maskedZeroCount);
    errors[ERR_ONE].resize(maskedCount);
    relErrors[ERR_ONE].resize(maskedCount);
    errorsDer[ERR_ONE].resize(maskedCount);

    int j = 0, k = 0;
    for (i = 0; i < levCount; i++)
    {
        curtime = levs[i].datetime;
        approxVector[cur]->GetLevels(curtime, 0, 1, &target, &filled, 0);

        if (filled == 1)
        {
            errors[ERR_ALL][i] = fabs(target - levs[i].level);
            relErrors[ERR_ALL][i] = errors[ERR_ALL][i] / levs[i].level;
        }
        else
            errors[ERR_ALL][i] = 0;

        approxVector[cur]->GetLevels(curtime, 0, 1, &target, &filled, 1);

        if (filled == 1)
            errorsDer[ERR_ALL][i] = fabs(target - refDerivations[i]);
        else
            errorsDer[ERR_ALL][i] = 0;

        if ((1 << (7 - i % 8)) & appCurrentTestMask)
        {
            errors[ERR_ONE][j] = errors[ERR_ALL][i];
            relErrors[ERR_ONE][j] = relErrors[ERR_ALL][i];
            errorsDer[ERR_ONE][j] = errorsDer[ERR_ALL][i];
            j++;
        }
        else
        {
            errors[ERR_ZERO][k] = errors[ERR_ALL][i];
            relErrors[ERR_ZERO][k] = relErrors[ERR_ALL][i];
            errorsDer[ERR_ZERO][k] = errorsDer[ERR_ALL][i];
            k++;
        }
    }

    floattype minval, maxval, med, fquart, tquart, avg, deviation;

    CalculateStats(errors[ERR_ALL], avg, minval, fquart, med, tquart, maxval, deviation);
    std::cout << avg << ";" << minval << ";" << fquart << ";" << med << ";" << tquart << ";" << maxval << ";" << deviation;
    std::cout << ";";
    CalculateStats(errors[ERR_ZERO], avg, minval, fquart, med, tquart, maxval, deviation);
    std::cout << avg << ";" << minval << ";" << fquart << ";" << med << ";" << tquart << ";" << maxval << ";" << deviation;
    std::cout << ";";
    CalculateStats(errors[ERR_ONE], avg, minval, fquart, med, tquart, maxval, deviation);
    std::cout << avg << ";" << minval << ";" << fquart << ";" << med << ";" << tquart << ";" << maxval << ";" << deviation;
    std::cout << ";";

    //rel
    CalculateStats(relErrors[ERR_ALL], avg, minval, fquart, med, tquart, maxval, deviation);
    std::cout << avg << ";" << minval << ";" << fquart << ";" << med << ";" << tquart << ";" << maxval << ";" << deviation;
    std::cout << ";";
    CalculateStats(relErrors[ERR_ZERO], avg, minval, fquart, med, tquart, maxval, deviation);
    std::cout << avg << ";" << minval << ";" << fquart << ";" << med << ";" << tquart << ";" << maxval << ";" << deviation;
    std::cout << ";";
    CalculateStats(relErrors[ERR_ONE], avg, minval, fquart, med, tquart, maxval, deviation);
    std::cout << avg << ";" << minval << ";" << fquart << ";" << med << ";" << tquart << ";" << maxval << ";" << deviation;
}

void CalculateAndPrintStats(std::vector<CGlucoseLevels*> &levelsVector, std::vector<CCommonApprox*> &approxVector)
{
    std::cout << "mask;segment;";
    std::cout << "avg(all);min(all);1Q(all);med(all);3Q(all);max(all);dev(all)";
    std::cout << ";";
    std::cout << "avg(0);min(0);1Q(0);med(0);3Q(0);max(0);dev(0)";
    std::cout << ";";
    std::cout << "avg(1);min(1);1Q(1);med(1);3Q(1);max(1);dev(1)";
    std::cout << ";";
    std::cout << "ravg(all);rmin(all);r1Q(all);rmed(all);r3Q(all);rmax(all);rdev(all)";
    std::cout << ";";
    std::cout << "ravg(0);rmin(0);r1Q(0);rmed(0);r3Q(0);rmax(0);rdev(0)";
    std::cout << ";";
    std::cout << "ravg(1);rmin(1);r1Q(1);rmed(1);r3Q(1);rmax(1);rdev(1)";

    for (size_t segment = 0; segment < levelsVector.size(); segment++)
    {
        for (size_t mask = 1; mask < 0x100; mask++)
        {
            appCurrentTestMask = (uint8_t)mask;
            std::cout << std::endl << (int)appCurrentTestMask << ";" << (int)segment << ";";
            CalculateAndPrintStatsForMask(segment, levelsVector, approxVector);
        }
    }

    std::cout << std::endl;

    if (appApproxMethod == apxmQuadraticSpline || appApproxMethod == apxmAkimaSpline || appApproxMethod == apxmCatmullRomSpline || appApproxMethod == apxmCubicHermiteSpline)
        std::cout << "The curve is continuously differentiable (C1 continuous)";
    else // no such case in this semestral work
        std::cout << "The curve is NOT continuously differentiable (not C1 continuous)";
}
