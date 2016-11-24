#include "../approx/src/GlucoseLevels.h"
#include "../loaders/SQLiteLoader.h"
#include "../loaders/Conversion.h"

#include <iostream>
#include <string>
#include <list>

SQLiteLoader::SQLiteLoader()
{
    db = nullptr;
}

HRESULT SQLiteLoader::InitLoader(LoaderParams params)
{
    int res;

    if (params.type != LoaderType::SQLiteLoaderType)
        return S_FALSE;

    // test if file exists
    FILE* test;
    fopen_s(&test, params.params.SQLiteLoaderParams.filename, "r");
    if (!test)
        return S_FALSE;
    fclose(test);

    res = sqlite3_open(params.params.SQLiteLoaderParams.filename, &db);
    if (res != 0)
    {
        sqlite3_close(db);
        return S_FALSE;
    }

    return S_OK;
}

HRESULT SQLiteLoader::LoadGlucoseLevels(std::vector<CGlucoseLevels*> &target)
{
    int res;
    int rows, columns;
    char* error;
    char **results = nullptr;
    int i;
    TGlucoseLevel* level;

    if (!db)
        return S_FALSE;

    std::list<int> segmentIds;

    // select all available segments with values
    std::string sql = "SELECT DISTINCT(segmentid) FROM measuredvalue;";
    res = sqlite3_get_table(db, sql.c_str(), &results, &rows, &columns, &error);
    if (res != 0)
    {
        std::cerr << "Error executing SQLite3 query: " << sqlite3_errmsg(db) << std::endl;
        sqlite3_free(error);
        return S_FALSE;
    }

    // fetch segment IDs to list for further use
    for (i = 1; i <= rows; ++i)
        segmentIds.push_back(std::stoi(results[i * columns + 0]));
    sqlite3_free_table(results);

    // select all values from each segment, and construct CGlucoseLevels instance
    for (int &segmentId : segmentIds)
    {
        // we acquire julianday, convert it to unix timestamp in milliseconds for further conversion to rational time

        sql = std::string("SELECT (JULIANDAY(measuredat) - 2440587.5) * 86400000.0 AS time, ist FROM measuredvalue WHERE ist IS NOT NULL AND segmentid = ") + std::to_string(segmentId) + " ORDER BY time ASC;";
        res = sqlite3_get_table(db, sql.c_str(), &results, &rows, &columns, &error);
        if (res != 0)
        {
            std::cerr << "Error executing SQLite3 query: " << sqlite3_errmsg(db) << std::endl;
            sqlite3_free(error);
            continue;
        }

        // create instance, set level count and retrieve pointer to array (vector)
        CGlucoseLevels* levs = new CGlucoseLevels();
        levs->SetLevelsCount(rows);
        levs->GetLevels(&level);

        // start from 1 due to sqlite3 including header in result
        for (i = 1; i <= rows; ++i)
        {
            level[i - 1].datetime = QDateTime2RatTime((int64_t)std::stod(results[i * columns + 0]));
            level[i - 1].level = std::stod(results[i * columns + 1]);
        }

        sqlite3_free_table(results);

        target.push_back(levs);
    }

    return S_OK;
}

HRESULT SQLiteLoader::Finalize()
{
    if (db)
        sqlite3_close(db);

    return S_OK;
}
