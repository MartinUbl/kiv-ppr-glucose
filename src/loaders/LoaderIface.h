#pragma once

#include <vector>

#include "../common/rtl/hresult.h"

class CGlucoseLevels;

// implemented loader types
enum class LoaderType
{
    SQLiteLoaderType,   // SQLite3 loader
};

// parameters passed to loader
struct LoaderParams
{
    // type of loader, see enum LoaderType
    LoaderType type;
    // parameters union, may contain different types of params depending on loader type
    union
    {
        struct
        {
            const char* filename;
        } SQLiteLoaderParams;
    } params;
};

// Loader interface class
class ILoader
{
    public:
        // initializes loader using supplied parameters
        virtual HRESULT InitLoader(LoaderParams params) = 0;
        // retrieves glucose levels and parses them into vector
        virtual HRESULT LoadGlucoseLevels(std::vector<CGlucoseLevels*> &target) = 0;
        // finalizes loader cycle (cleanup)
        virtual HRESULT Finalize() = 0;
};
