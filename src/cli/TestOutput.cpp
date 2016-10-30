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

        for (size_t i = 0; i < 48; i++)
        {
            if (consoleout)
                std::cout << "Original time: " << levs[i].datetime - base << ", value: " << levs[i].level << std::endl;

            ofs << "<circle stroke=\"#FF0000\" fill=\"#FFFFFF\" r=\"5\" cx=\"" << (levs[i].datetime) * stretchCoefX - base << "\" cy=\"" << (invertBase - levs[i].level) * stretchCoefY << "\" />" << std::endl;
        }

        floattype tim = starttime;
        for (size_t i = 0; i < count; i++)
        {
            if (consoleout)
                std::cout << "Time: " << tim - base << ", value: " << levels[i] << std::endl;

            ofs << "<circle stroke=\"#0000FF\" fill=\"#FFFFFF\" r=\"2\" cx=\"" << (tim) * stretchCoefX - base << "\" cy=\"" << (invertBase - levels[i]) * stretchCoefY << "\" />" << std::endl;
            tim += step;
        }

        ofs << "</svg>";
    }
}
