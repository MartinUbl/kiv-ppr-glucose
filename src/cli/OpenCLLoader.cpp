#include "OpenCLLoader.h"
#include "../common/iface/ApproxIface.h"
#include "appconfig.h"

#include <iostream>
#include <map>

#define MAX_CL_SOURCE_SIZE (0x100000)

std::map<size_t, clProgramRecord> clProgramMap;

bool _clSupportsDouble = false;

bool clSupportsDouble()
{
    return _clSupportsDouble;
}

bool LoadCLPrograms()
{
    FILE* fp;
    char *source_str;
    size_t source_size;
    cl_int ret;

    // for each program listed in config
    for (int i = 0; i < sizeof(clPrograms) / sizeof(clProgramLocator); i++)
    {
        // attempt to open file from current workdir
        fopen_s(&fp, clPrograms[i].clFile, "r");
        if (!fp)
        {
            // if no such file, try debug file path
            std::string dbgPath = OPENCL_DEBUG_LOADER_BASE;
            dbgPath += clPrograms[i].clFile;
            fopen_s(&fp, dbgPath.c_str(), "r");

            if (!fp)
            {
                std::cerr << "Failed to load kernel code from " << clPrograms[i].clFile << std::endl;
                return false;
            }
        }

        // load program
        source_str = new char[MAX_CL_SOURCE_SIZE];
        source_size = fread(source_str, 1, MAX_CL_SOURCE_SIZE, fp);
        fclose(fp);

        // create context, initialize command queue
        cl_platform_id platform_id = NULL;
        cl_device_id device_id = NULL;
        cl_uint ret_num_devices;
        cl_uint ret_num_platforms;
        ret = clGetPlatformIDs(1, &platform_id, &ret_num_platforms);
        ret = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_DEFAULT, 1, &device_id, &ret_num_devices);
        clProgramMap[clPrograms[i].apxMethod].context = clCreateContext(NULL, 1, &device_id, NULL, NULL, &ret);

        clProgramMap[clPrograms[i].apxMethod].commandQueue = clCreateCommandQueue(clProgramMap[clPrograms[i].apxMethod].context, device_id, 0, &ret);

        // create and build program
        clProgramMap[clPrograms[i].apxMethod].prog = clCreateProgramWithSource(clProgramMap[clPrograms[i].apxMethod].context, 1, (const char **)&source_str, (const size_t *)&source_size, &ret);
        ret = clBuildProgram(clProgramMap[clPrograms[i].apxMethod].prog, 1, &device_id, NULL, NULL, NULL);

        // determine double precision availability
        cl_device_fp_config double_cfg;
        clGetDeviceInfo(device_id, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof(double_cfg), &double_cfg, NULL);
        _clSupportsDouble = (double_cfg != 0);

        // cleanup
        delete[] source_str;
    }

    return true;
}

void FinalizeCLPrograms()
{
    // for each program loaded
    for (std::map<size_t, clProgramRecord>::iterator itr = clProgramMap.begin(); itr != clProgramMap.end(); ++itr)
    {
        // flush command queue
        clFlush(itr->second.commandQueue);
        clFinish(itr->second.commandQueue);

        // release resources from structure
        clReleaseProgram(itr->second.prog);
        clReleaseCommandQueue(itr->second.commandQueue);
        clReleaseContext(itr->second.context);
    }

    clProgramMap.clear();
}

clProgramRecord* GetCLProgramRecord(size_t apxMethod)
{
    if (clProgramMap.find(apxMethod) == clProgramMap.end())
        return nullptr;

    return &clProgramMap[apxMethod];
}
