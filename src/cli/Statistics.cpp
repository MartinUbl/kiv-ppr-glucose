#include "Statistics.h"
#include "appconfig.h"
#include <iostream>

void CalculateAndPrintStats(std::vector<CGlucoseLevels*> &levelsVector, std::vector<CCommonApprox*> &approxVector)
{
    size_t cur = 0;
    appCurrentTestMask = 22;
    TGlucoseLevel* levs;
    size_t levCount;
    size_t filled;
    floattype target, curtime;

    levelsVector[cur]->GetLevels(&levs);
    levelsVector[cur]->GetLevelsCount(&levCount);

    for (int i = 0; i < levCount; i++)
    {
        curtime = levs[i].datetime;
        approxVector[cur]->GetLevels(curtime, 0, 1, &target, &filled, 0);

        if (filled == 1)
            std::cout << "Err: " << (target - levs[i].level) << std::endl;
    }

    std::cout << std::endl;
}
