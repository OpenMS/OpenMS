// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzMLSqliteSwathHandler.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/FORMAT/SqliteConnector.h>
#include <sqlite3.h>

namespace OpenMS::Internal
{


    namespace Sql = Internal::SqliteHelper;

    std::vector<OpenSwath::SwathMap> MzMLSqliteSwathHandler::readSwathWindows()
    {
      std::vector<OpenSwath::SwathMap> swath_maps;
      SqliteConnector conn(filename_);
      sqlite3_stmt * stmt;

      std::string select_sql;
      select_sql = "SELECT " \
                    "DISTINCT(ISOLATION_TARGET)," \
                    "ISOLATION_TARGET - ISOLATION_LOWER," \
                    "ISOLATION_TARGET + ISOLATION_UPPER " \
                    "FROM PRECURSOR " \
                    "INNER JOIN SPECTRUM ON SPECTRUM_ID = SPECTRUM.ID " \
                    "WHERE MSLEVEL == 2 "\
                    ";";

      conn.prepareStatement(&stmt, select_sql);
      sqlite3_step( stmt );

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        OpenSwath::SwathMap map;
        Sql::extractValue<double>(&map.center, stmt, 0);
        Sql::extractValue<double>(&map.lower, stmt, 1);
        Sql::extractValue<double>(&map.upper, stmt, 2);
        swath_maps.push_back(map);
        sqlite3_step( stmt );
      }

      // free memory
      sqlite3_finalize(stmt);

      return swath_maps;
    }

    std::vector<int> MzMLSqliteSwathHandler::readMS1Spectra()
    {
      std::vector< int > indices;
      SqliteConnector conn(filename_);
      sqlite3_stmt * stmt;

      std::string select_sql;
      select_sql = "SELECT ID " \
                   "FROM SPECTRUM " \
                   "WHERE MSLEVEL == 1;";

      conn.prepareStatement(&stmt, select_sql);
      sqlite3_step(stmt);

      while (sqlite3_column_type(stmt, 0) != SQLITE_NULL)
      {
        indices.push_back(sqlite3_column_int(stmt, 0));
        sqlite3_step(stmt);
      }

      // free memory
      sqlite3_finalize(stmt);

      return indices;
    }

    std::vector<int> MzMLSqliteSwathHandler::readSpectraForWindow(const OpenSwath::SwathMap& swath_map)
    {
      std::vector< int > indices;
      const double center = swath_map.center;

      SqliteConnector conn(filename_);
      sqlite3_stmt * stmt;

      String select_sql = "SELECT " \
                          "SPECTRUM_ID " \
                          "FROM PRECURSOR " \
                          "WHERE ISOLATION_TARGET BETWEEN ";

      select_sql += String(center - 0.01) + " AND " + String(center + 0.01) + ";";
      conn.prepareStatement(&stmt, select_sql);
      sqlite3_step(stmt);

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        indices.push_back(sqlite3_column_int(stmt, 0));
        sqlite3_step(stmt);
      }

      // free memory
      sqlite3_finalize(stmt);

      return indices;
    }

} // namespace OpenMS  // namespace Internal

