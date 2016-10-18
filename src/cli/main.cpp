#include <iostream>
#include <iomanip>
#include <vector>
#include "../approx/src/GlucoseLevels.h"
#include "../loaders/SQLiteLoader.h"
#include "../approx/src/ApproxQuadraticSpline.h"
#include "../approx/src/ApproxAkimaSpline.h"
#include "appconfig.h"

#include "TestOutput.h"

// concurrency used for calculations
ConcurrencyType appConcurrency = ConcurrencyType::ct_serial;
// selected approximation method
size_t appApproxMethod = 0;
// worker count (if the method needs it)
size_t appWorkerCount = 1;

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

    // TODO: load approx method from command line parameters
    appApproxMethod = apxmAkimaSpline;

    // TODO: load concurrency type from command line parameters
    appConcurrency = ConcurrencyType::ct_parallel_threads;

    // TODO: load worker count from command line parameters
    appWorkerCount = 4;

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

    const floattype timestart = 0.0;

    clock_t tmStart = clock();

    std::cout << "Processing " << vec.size() << " segments..." << std::endl;
    for (size_t i = 0; i < vec.size(); i++)
    {
        std::cout << "Processing segment " << (i+1) << " / " << vec.size() << "..." << std::endl;

        // reduce times to not lose precision so rapidly
        floattype reducedBy = reduceLevels(vec[i]);

        CCommonApprox* apx;
        // approximate!
        if (appApproxMethod == apxmQuadraticSpline)
            apx = new ApproxQuadraticSpline(vec[i]);
        else if (appApproxMethod == apxmAkimaSpline)
            apx = new ApproxAkimaSpline(vec[i]);
        else
            break;

        apx->Approximate(nullptr);

        // some testing values for output
        size_t levcount = 768;
        size_t filled;
        floattype* levels = new floattype[levcount];

        const floattype granularity = 16.0;
        const floattype step = 5.0 / (24.0*60.0*granularity);

        // retrieve approximated levels
        apx->GetLevels(timestart, step, levcount, levels, &filled, 0);

        /*
        if (i == 0)
        {
            TGlucoseLevel* levs;
            vec[i]->GetLevels(&levs);

            // visualize to some easy format, SVG should be nice
            VisualizeSVG("test.svg", timestart, step, levcount, levels, levs, false);
        }
        */

        delete apx;
    }

    clock_t tmTotal = (clock_t)(1000.0f * (float)(clock() - tmStart) / (float)CLOCKS_PER_SEC);

    std::cout << "Done. Elapsed: " << tmTotal << "ms" << std::endl;

    ldr.Finalize();

    return 0;
}
