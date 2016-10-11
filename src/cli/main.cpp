#include <iostream>
#include <iomanip>
#include <vector>
#include "../approx/src/GlucoseLevels.h"
#include "../loaders/SQLiteLoader.h"
#include "../approx/src/ApproxQuadraticSpline.h"
#include "../approx/src/ApproxAkimaSpline.h"

#include "TestOutput.h"

// reduces all times by minimum of all times - causes less precision loss
static floattype reduceLevels(CGlucoseLevels* lvls)
{
    TGlucoseLevel* lev;
    size_t cnt;
    lvls->GetLevels(&lev);
    lvls->GetLevelsCount(&cnt);

    floattype reducer = lev[0].datetime;

    for (size_t i = 0; i < cnt; i++)
        lev[i].datetime -= reducer;

    return reducer;
}

int main(int argc, char** argv)
{
    // TODO

    // Load parameters - fow now use just SQLite loader
    SQLiteLoader ldr;
    HRESULT res = ldr.InitLoader({
        LoaderType::SQLiteLoaderType,
        {
            "../../data/direcnet.sqlite"
        }
    });

    std::vector<CGlucoseLevels*> vec;
    res = ldr.LoadGlucoseLevels(vec);

    // reduce times to not lose precision so rapidly
    floattype reducedBy = reduceLevels(vec[0]);
    const floattype timestart = 0.0;

    std::cout << std::setprecision(20);

    // approximate!
    //ApproxQuadraticSpline apx(vec[0]);
    ApproxAkimaSpline apx(vec[0]);
    apx.Approximate(nullptr);

    // some testing values for output
    size_t levcount = 768;
    size_t filled;
    floattype* levels = new floattype[levcount];

    const floattype granularity = 16.0;
    const floattype step = 5.0 / (24.0*60.0*granularity);

    // retrieve approximated levels
    apx.GetLevels(timestart, step, levcount, levels, &filled, 0);

    TGlucoseLevel* levs;
    vec[0]->GetLevels(&levs);

    // visualize to some easy format, SVG should be nice
    VisualizeSVG("test.svg", timestart, step, levcount, levels, levs, true);

    ldr.Finalize();

    return 0;
}
