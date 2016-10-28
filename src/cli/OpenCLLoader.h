#pragma once

#ifndef CL_USE_DEPRECATED_OPENCL_2_0_APIS
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#endif
#ifndef CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#endif
#include <CL/cl.h>

#define OPENCL_DEBUG_LOADER_BASE "../../src/approx/src/"

// OpenCL prepared program structure
typedef struct
{
    // OpenCL program
    cl_program prog;
    // program context
    cl_context context;
    // program command queue
    cl_command_queue commandQueue;
} clProgramRecord;

// does local device support double precision?
bool clSupportsDouble();

// loads all available OpenCL programs
bool LoadCLPrograms();
// finalizes loaded OpenCL programs
void FinalizeCLPrograms();
// retrieves OpenCL program record (if available)
clProgramRecord* GetCLProgramRecord(size_t apxMethod);
