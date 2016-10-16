#pragma once

enum class ConcurrencyType
{
    ct_serial,                  // no parallelization at all
    cr_parallel_threads,        // use C++ std::thread
    ct_parallel_amp,            // use C++AMP
};

// concurrency used for calculations
static ConcurrencyType appConcurrency = ConcurrencyType::ct_serial;
// selected approximation method
static size_t appApproxMethod = 0;
// worker count (if the method needs it)
static size_t appWorkerCount = 1;
