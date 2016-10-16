#pragma once

enum class ConcurrencyType
{
    ct_serial,                  // no parallelization at all
    ct_parallel_threads,        // use C++ std::thread
    ct_parallel_amp,            // use C++AMP
};

// concurrency used for calculations
extern ConcurrencyType appConcurrency;
// selected approximation method
extern size_t appApproxMethod;
// worker count (if the method needs it)
extern size_t appWorkerCount;
