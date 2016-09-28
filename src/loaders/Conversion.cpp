#include "Conversion.h"

double QDateTime2RatTime(const int64_t unixepochtime)
{
    const int64_t diffFrom1970To1900 = 2209161600000;
    const double MSecsPerDay = 24.0*60.0*60.0*1000.0;
    const double InvMSecsPerDay = 1.0 / MSecsPerDay;

    int64_t diff = unixepochtime + diffFrom1970To1900;
    return ((double)diff)*InvMSecsPerDay;
}
