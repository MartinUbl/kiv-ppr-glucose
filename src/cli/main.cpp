#include <iostream>
#include <iomanip>
#include <vector>
#include "../approx/src/GlucoseLevels.h"
#include "../loaders/SQLiteLoader.h"
#include "../approx/src/ApproxQuadraticSpline.h"
#include "../approx/src/ApproxAkimaSpline.h"
#include "appconfig.h"
#include "OpenCLLoader.h"

#include "TestOutput.h"

// concurrency used for calculations
ConcurrencyType appConcurrency = ConcurrencyType::ct_serial;
// selected approximation method
size_t appApproxMethod = 0;
// worker count (if the method needs it)
size_t appWorkerCount = 1;
// console output enabled?
bool appSilentMode = false;

// should we create testing SVG output? (1 = true, anything else = false)
#define DEBUG_SVG_PRINT 0

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
    // TODO: load approx method from command line parameters
    appApproxMethod = apxmQuadraticSpline;

    // TODO: load concurrency type from command line parameters
    appConcurrency = ConcurrencyType::ct_parallel_opencl;

    // TODO: load worker count from command line parameters
    appWorkerCount = 4;

    // TODO: load silent mode parameter from command line parameters
    appSilentMode = false;

    // Load parameters - fow now use just SQLite loader
    SQLiteLoader ldr;
    HRESULT res = ldr.InitLoader({
        LoaderType::SQLiteLoaderType,
        {
            "../../data/direcnet.sqlite"
        }
    });

    if (!appSilentMode)
        std::cout << "Loading values from storage..." << std::endl;

    std::vector<CGlucoseLevels*> vec;
    res = ldr.LoadGlucoseLevels(vec);

    const floattype timestart = 0.0;

    // if we use OpenCL concurrency type, preload programs before actual calculation
    if (appConcurrency == ConcurrencyType::ct_parallel_opencl)
    {
        if (!appSilentMode)
            std::cout << "Loading OpenCL kernel files..." << std::endl;

        if (!LoadCLPrograms())
        {
            std::cerr << "Failed to load OpenCL kernel files, exiting." << std::endl;
            return 1;
        }
    }

    CCommonApprox* apx = new ApproxAkimaSpline(vec[0]);
    apx->Approximate(nullptr);

    clock_t tmStart = clock();

    uint64_t totalCount = 0;

    if (!appSilentMode)
        std::cout << "Processing " << vec.size() << " segments..." << std::endl;

    for (size_t i = 0; i < vec.size(); i++)
    {
        if (!appSilentMode)
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

#if DEBUG_SVG_PRINT == 1
        // some testing values for output
        size_t levcount = 4768;
        size_t filled;
        floattype* levels = new floattype[levcount];

        const floattype granularity = 16.0;
        const floattype step = 5.0 / (24.0*60.0*granularity);

        // retrieve approximated levels
        apx->GetLevels(timestart, step, levcount, levels, &filled, 0);

        if (i == 0)
        {
            TGlucoseLevel* levs;
            vec[i]->GetLevels(&levs);

            // visualize to some easy format, SVG should be nice
            VisualizeSVG("test.svg", timestart, step, levcount, levels, levs, false);
        }
#endif

        delete apx;
    }

    clock_t tmTotal = (clock_t)(1000.0f * (float)(clock() - tmStart) / (float)CLOCKS_PER_SEC);

    if (!appSilentMode)
        std::cout << "Done. Elapsed: " << tmTotal << "ms" << std::endl;

    ldr.Finalize();

    // if we use OpenCL concurrency type, erase preloaded programs
    if (appConcurrency == ConcurrencyType::ct_parallel_opencl)
    {
        if (!appSilentMode)
            std::cout << "Erasing OpenCL kernel files..." << std::endl;

        FinalizeCLPrograms();
    }

    return 0;
}
