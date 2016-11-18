#include <iostream>
#include <fstream>
#include <iomanip>
#include "TestOutput.h"

void VisualizeSVG(const char* filename, floattype starttime, floattype step, size_t count, floattype* levels, TGlucoseLevel* levs, bool consoleout)
{
    const floattype stretchCoefX = 10000.0;
    const floattype stretchCoefY = 100.0;
    const floattype invertBase = 12.0;

    // fstream scope
    {
        std::ofstream ofs(filename);

        ofs << "<svg xmlns=\"http://www.w3.org/2000/svg\" stroke-width=\"1\">" << std::endl;
        ofs << std::setprecision(20);
        if (consoleout)
            std::cout << std::setprecision(20);

        floattype base = levs[0].datetime * stretchCoefX;

        floattype tim = starttime;
        floattype px, py;
        for (size_t i = 0; i < count; i++)
        {
            if (consoleout)
                std::cout << "Time: " << tim - base << ", value: " << levels[i] << std::endl;

            // skip the first, just store its value later to be able to draw line in next iteration
            // also skip too large jumps
            if (i > 0 && fabs(py - (invertBase - levels[i]) * stretchCoefY) < 8.0 * stretchCoefY)
                ofs << "<line stroke=\"#0000FF\" stroke-width=\"2\" x1=\"" << px << "\" y1=\"" << py << "\" x2=\"" << (tim) * stretchCoefX - base << "\" y2=\"" << (invertBase - levels[i]) * stretchCoefY << "\" />" << std::endl;

            px = (tim)* stretchCoefX - base;
            py = (invertBase - levels[i]) * stretchCoefY;
            tim += step;
        }

        for (size_t i = 0; i < 48; i++)
        {
            if (consoleout)
                std::cout << "Original time: " << levs[i].datetime - base << ", value: " << levs[i].level << std::endl;

            ofs << "<circle stroke=\"#FF0000\" fill=\"#FF0000\" r=\"3\" cx=\"" << (levs[i].datetime) * stretchCoefX - base << "\" cy=\"" << (invertBase - levs[i].level) * stretchCoefY << "\" />" << std::endl;
        }

        ofs << "</svg>";
    }
}
