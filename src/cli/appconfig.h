#pragma once

enum class ConcurrencyType
{
    ct_serial,                  // no parallelization at all
    ct_parallel_threads,        // use C++ std::thread
    ct_parallel_amp_gpu,        // use C++AMP on GPU
    ct_parallel_tbb,            // use Intel Threading Building Blocks (TBB)
    ct_parallel_opencl,         // use OpenCL
};

// concurrency used for calculations
extern ConcurrencyType appConcurrency;
// selected approximation method
extern size_t appApproxMethod;
// worker count (if the method needs it)
extern size_t appWorkerCount;
// console output enabled?
extern bool appSilentMode;
