#pragma once

#include "../loaders/LoaderIface.h"

#include "../../dep/sqlite/sqlite3.h"
#pragma comment(lib, "../../dep/sqlite/sqlite3.lib")

class SQLiteLoader : public ILoader
{
    public:
        SQLiteLoader();

        HRESULT InitLoader(LoaderParams params);
        HRESULT LoadGlucoseLevels(std::vector<CGlucoseLevels*> &target);
        HRESULT Finalize();

    protected:
        sqlite3 *db;
};
