#pragma once

#include <string>

// maximum worker count (for threaded version)
#define MAX_APP_WORKER_COUNT 32

enum class ConcurrencyType
{
    ct_none,                    // dummy value used just for initialization
    ct_serial,                  // no parallelization at all
    ct_parallel_threads,        // use C++ std::thread
    ct_parallel_amp_gpu,        // use C++AMP on GPU
    ct_parallel_tbb,            // use Intel Threading Building Blocks (TBB)
    ct_parallel_opencl,         // use OpenCL
};

enum class LoaderType;

// concurrency used for calculations
extern ConcurrencyType appConcurrency;
// selected approximation method
extern size_t appApproxMethod;
// worker count (if the method needs it)
extern size_t appWorkerCount;
// console output enabled?
extern bool appSilentMode;
// input filename
extern std::string appInputFilename;
// loader type
extern LoaderType appLoaderType;

// current testing mask
extern uint8_t appCurrentTestMask;

typedef struct
{
    size_t apxMethod;
    const char* clFile;
} clProgramLocator;

// available OpenCL kernel files
const clProgramLocator clPrograms[] = {
    { apxmQuadraticSpline, "ApproxQuadraticSpline.cl" },
    { apxmAkimaSpline, "ApproxAkimaSpline.cl" }
};
