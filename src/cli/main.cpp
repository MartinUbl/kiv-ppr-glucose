#include <iostream>
#include <iomanip>
#include <vector>
#include "../approx/src/GlucoseLevels.h"
#include "../loaders/SQLiteLoader.h"
#include "../approx/src/ApproxQuadraticSpline.h"
#include "../approx/src/ApproxAkimaSpline.h"
#include "appconfig.h"
#include "OpenCLLoader.h"
#include "Statistics.h"

#include "TestOutput.h"

// concurrency used for calculations
ConcurrencyType appConcurrency;
// selected approximation method
size_t appApproxMethod;
// worker count (if the method needs it)
size_t appWorkerCount;
// console output enabled?
bool appSilentMode;
// input filename
std::string appInputFilename;
// loader type
LoaderType appLoaderType;

// current testing mask
uint8_t appCurrentTestMask;

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

static void printHelp(const char* progname)
{
    // extract program name
    std::string fname = progname;
    size_t p = fname.find_last_of('\\');
    if (p != std::string::npos)
        fname = fname.substr(p + 1);
    p = fname.find_last_of('/');
    if (p != std::string::npos)
        fname = fname.substr(p + 1);

    std::cout << "Usage: " << fname.c_str() << " -i <input filename> -m <method name> -c <concurrency type> [-id <storage driver>] [-w <worker count>] [-s] " << std::endl;
    std::cout << std::endl;
    std::cout << "Available methods: quadratic (q) - quadratic spline" << std::endl;
    std::cout << "                   akima (a) - akima spline" << std::endl;
    std::cout << "                   ??? (?) - ?" << std::endl;
    std::cout << std::endl;
    std::cout << "Available concurrency: serial (s) - serial version" << std::endl;
    std::cout << "                       threads (t) - standard C++ threads" << std::endl;
    std::cout << "                       tbb - Intel TBB" << std::endl;
    std::cout << "                       opencl (cl) - OpenCL" << std::endl;
    std::cout << std::endl;
    std::cout << "Available input drivers: sqlite - SQLite3 driver" << std::endl;
    std::cout << std::endl;
}

size_t parseApproxMethod(const char* str)
{
    if (strcmp(str, "quadratic") == 0 || strcmp(str, "q") == 0)
        return apxmQuadraticSpline;
    else if (strcmp(str, "akima") == 0 || strcmp(str, "a") == 0)
        return apxmAkimaSpline;

    return 0;
}

ConcurrencyType parseConcurrencyType(const char* str)
{
    if (strcmp(str, "serial") == 0 || strcmp(str, "s") == 0)
        return ConcurrencyType::ct_serial;
    else if (strcmp(str, "tbb") == 0)
        return ConcurrencyType::ct_parallel_tbb;
    else if (strcmp(str, "threads") == 0 || strcmp(str, "t") == 0)
        return ConcurrencyType::ct_parallel_threads;
    else if (strcmp(str, "opencl") == 0 || strcmp(str, "cl") == 0)
        return ConcurrencyType::ct_parallel_opencl;
    //else if (strcmp(str, "cppamp") == 0 || strcmp(str, "amp") == 0)
    //    return ConcurrencyType::ct_parallel_amp_gpu;

    return ConcurrencyType::ct_none;
}

int parseCLIArgs(int argc, char** argv)
{
    // 1 program name, 3 parameters, 3 values
    if (argc < 7)
        return 1;

    // load defaults
    appInputFilename = "";
    appApproxMethod = 0;
    appConcurrency = ConcurrencyType::ct_none;
    appLoaderType = LoaderType::SQLiteLoaderType;
    appSilentMode = false;
    appWorkerCount = 1;

    for (int i = 1; i < argc; i++)
    {
        if (strcmp(argv[i], "-i") == 0)
        {
            if (i < argc - 1)
            {
                appInputFilename = argv[i + 1];
                i++;
            }
            else
                std::cerr << "No filename specified after -i parameter." << std::endl;
        }
        else if (strcmp(argv[i], "-m") == 0)
        {
            if (i < argc - 1)
            {
                appApproxMethod = parseApproxMethod(argv[i + 1]);
                if (appApproxMethod == 0)
                    std::cerr << "Invalid approximation/interpolation method '" << argv[i + 1] << "' specified." << std::endl;
            }
            else
                std::cerr << "No approximation/interpolation method specified after -m parameter." << std::endl;
        }
        else if (strcmp(argv[i], "-c") == 0)
        {
            if (i < argc - 1)
            {
                appConcurrency = parseConcurrencyType(argv[i + 1]);
                if (appConcurrency == ConcurrencyType::ct_none)
                    std::cerr << "Invalid concurrency '" << argv[i + 1] << "' specified." << std::endl;
            }
            else
                std::cerr << "No concurrency type specified after -c parameter." << std::endl;
        }
        else if (strcmp(argv[i], "-id") == 0)
        {
            if (i < argc - 1)
            {
                appLoaderType = ((strcmp(argv[i + 1], "sqlite") == 0) ? LoaderType::SQLiteLoaderType : LoaderType::UnknownLoaderType);
                if (appLoaderType == LoaderType::UnknownLoaderType)
                    std::cerr << "Invalid storage driver '" << argv[i + 1] << "' specified." << std::endl;
            }
            else
                std::cerr << "No storage driver specified after -id parameter." << std::endl;
        }
        else if (strcmp(argv[i], "-w") == 0)
        {
            if (i < argc - 1)
            {
                char* endptr;
                appWorkerCount = std::strtoul(argv[i + 1], &endptr, 10);
                if (appWorkerCount <= 0 || appWorkerCount > MAX_APP_WORKER_COUNT)
                    std::cerr << "Invalid worker count '" << argv[i + 1] << "' specified (must be number in range 1 - " << MAX_APP_WORKER_COUNT << "." << std::endl;
            }
            else
                std::cerr << "No worker count specified after -w parameter." << std::endl;
        }
        else if (strcmp(argv[i], "-s") == 0)
        {
            appSilentMode = true;
        }
    }

    if (appInputFilename == "" || appApproxMethod == 0 || appConcurrency == ConcurrencyType::ct_none)
    {
        std::cerr << "Not enough parameters specified" << std::endl;
        return 4;
    }

    return 0;
}

int main(int argc, char** argv)
{
    // parse input parameters
    if (parseCLIArgs(argc, argv) != 0)
    {
        printHelp(argv[0]);
        return 1;
    }

    // Load parameters - fow now use just SQLite loader
    ILoader* ldr = nullptr;

    if (appLoaderType == LoaderType::SQLiteLoaderType)
        ldr = new SQLiteLoader();
    else
    {
        std::cerr << "Unsupported loader type." << std::endl;
        return 2;
    }

    // instantiate loader
    HRESULT res = ldr->InitLoader({
        appLoaderType,
        {
            appInputFilename.c_str()
        }
    });

    // loader initialization succeeded?
    if (res == S_FALSE)
    {
        std::cerr << "Failed to load file " << appInputFilename.c_str() << " - invalid file format or file not found." << std::endl;
        return 3;
    }

    if (!appSilentMode)
    {
        std::cout << "Using approximator: ";
        if (appApproxMethod == apxmQuadraticSpline)
            std::cout << "quadratic spline";
        else if (appApproxMethod == apxmAkimaSpline)
            std::cout << "akima spline";
        else
            std::cout << "unknown";
        std::cout << std::endl;

        std::cout << "Using concurrency mode: ";
        if (appConcurrency == ConcurrencyType::ct_serial)
            std::cout << "none (serial)";
        else if (appConcurrency == ConcurrencyType::ct_parallel_threads)
            std::cout << "C++11 threads";
        else if (appConcurrency == ConcurrencyType::ct_parallel_opencl)
            std::cout << "OpenCL";
        else if (appConcurrency == ConcurrencyType::ct_parallel_tbb)
            std::cout << "Intel TBB";
        else if (appConcurrency == ConcurrencyType::ct_parallel_amp_gpu)
            std::cout << "C++AMP (GPU)";
        else
            std::cout << "unknown";
        std::cout << std::endl;
    }

    if (!appSilentMode)
        std::cout << "Loading values from storage..." << std::endl;

    // load values from storage
    std::vector<CGlucoseLevels*> vec;
    res = ldr->LoadGlucoseLevels(vec);

    size_t vecSize = vec.size();
    size_t i;

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

    if (!appSilentMode)
        std::cout << "Preparing approximation structures..." << std::endl;

    std::vector<CCommonApprox*> approxVect(vecSize);
    std::vector<floattype> reduceAmount(vecSize);
    for (i = 0; i < vecSize; i++)
    {
        // reduce times to not lose precision so rapidly
        reduceAmount[i] = reduceLevels(vec[i]);

        if (appApproxMethod == apxmQuadraticSpline)
            approxVect[i] = new ApproxQuadraticSpline(vec[i]);
        else if (appApproxMethod == apxmAkimaSpline)
            approxVect[i] = new ApproxAkimaSpline(vec[i]);
    }

    if (!appSilentMode)
        std::cout << "Processing " << vec.size() << " segments..." << std::endl;

    clock_t tmStart = clock();

    for (i = 0; i < vecSize; i++)
    {
        approxVect[i]->Approximate(nullptr);

#if DEBUG_SVG_PRINT == 1
        if (i == 0)
        {
            appCurrentTestMask = 154;
            // some testing values for output
            size_t levcount = 768;
            size_t filled;
            floattype* levels = new floattype[levcount];

            const floattype granularity = 16.0;
            const floattype step = 5.0 / (24.0*60.0*granularity);

            // retrieve approximated levels
            approxVect[i]->GetLevels(0.0, step, levcount, levels, &filled, 0);

            TGlucoseLevel* levs;
            vec[i]->GetLevels(&levs);

            // visualize to some easy format, SVG should be nice
            VisualizeSVG("test.svg", 0.0, step, levcount, levels, levs, false);
        }
#endif
    }

    clock_t tmTotal = (clock_t)(1000.0f * (float)(clock() - tmStart) / (float)CLOCKS_PER_SEC);

    if (!appSilentMode)
        std::cout << "Done. Elapsed: " << tmTotal << "ms" << std::endl;

    ldr->Finalize();

    // if we use OpenCL concurrency type, erase preloaded programs
    if (appConcurrency == ConcurrencyType::ct_parallel_opencl)
    {
        if (!appSilentMode)
            std::cout << "Erasing OpenCL kernel files..." << std::endl;

        FinalizeCLPrograms();
    }

    if (!appSilentMode)
        std::cout << std::endl << "Done calculating parameters. Calculating statistics..." << std::endl;

    CalculateAndPrintStats(vec, approxVect);

    return 0;
}
