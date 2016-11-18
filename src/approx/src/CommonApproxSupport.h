#pragma once

#include "../../common/iface/ApproxIface.h"

#ifndef M_PI
// Ludolph's number - PI
#define M_PI 3.14159265358979323846
#endif

/*
 * Solve cubic equation with only one root in range <0;1>
 */
floattype Cubic_IdentitySolve(floattype a, floattype b, floattype c, floattype d);
