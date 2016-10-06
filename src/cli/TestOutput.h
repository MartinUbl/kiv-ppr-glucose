#pragma once

#include "../common/iface/ApproxIface.h"

void VisualizeSVG(const char* filename, floattype starttime, floattype step, size_t count, floattype* levels, TGlucoseLevel* levs, bool consoleout = false);
