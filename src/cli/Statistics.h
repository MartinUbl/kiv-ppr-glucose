#pragma once

#include <vector>
#include "../approx/src/CommonApprox.h"
#include "../approx/src/GlucoseLevels.h"

void CalculateAndPrintStats(std::vector<CGlucoseLevels*> &levelsVector, std::vector<CCommonApprox*> &approxVector);
