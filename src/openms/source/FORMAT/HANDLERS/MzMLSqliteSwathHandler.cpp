// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

namespace OpenMS
{
  namespace Internal
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

  } // namespace Internal
} // namespace OpenMS

