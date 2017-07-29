// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <sqlite3.h>

namespace OpenMS
{
  namespace Internal
  {

    std::vector<OpenSwath::SwathMap> MzMLSqliteSwathHandler::readSwathWindows()
    {
      std::vector<OpenSwath::SwathMap> swath_maps;
      sqlite3 *db = openDB();
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

      sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      sqlite3_step( stmt );

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        OpenSwath::SwathMap map;

        if (sqlite3_column_type(stmt, 0) != SQLITE_NULL) map.center = sqlite3_column_double(stmt, 0);
        if (sqlite3_column_type(stmt, 1) != SQLITE_NULL) map.lower = sqlite3_column_double(stmt, 1);
        if (sqlite3_column_type(stmt, 2) != SQLITE_NULL) map.upper = sqlite3_column_double(stmt, 2);

        swath_maps.push_back(map);
        sqlite3_step( stmt );
      }

      // free memory and close connection
      sqlite3_finalize(stmt);
      sqlite3_close(db);

      return swath_maps;
    }

    std::vector<int> MzMLSqliteSwathHandler::readMS1Spectra()
    {
      std::vector< int > indices;
      sqlite3 *db = openDB();
      sqlite3_stmt * stmt;

      std::string select_sql;
      select_sql = "SELECT ID " \
                   "FROM SPECTRUM " \
                   "WHERE MSLEVEL == 1;";

      sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      sqlite3_step(stmt);

      while (sqlite3_column_type(stmt, 0) != SQLITE_NULL)
      {
        indices.push_back(sqlite3_column_int(stmt, 0));
        sqlite3_step(stmt);
      }

      // free memory and close connection
      sqlite3_finalize(stmt);
      sqlite3_close(db);

      return indices;
    }

    std::vector<int> MzMLSqliteSwathHandler::readSpectraForWindow(OpenSwath::SwathMap swath_map)
    {
      std::vector< int > indices;
      double center = swath_map.center;

      sqlite3 *db = openDB();
      sqlite3_stmt * stmt;

      std::string select_sql;
      select_sql = "SELECT " \
                    "SPECTRUM_ID " \
                    "FROM PRECURSOR " \
                    "WHERE ISOLATION_TARGET BETWEEN "; 

      select_sql += String(center - 0.01); 
      select_sql += " AND ";
      select_sql += String(center + 0.01); 
      select_sql += ";";

      sqlite3_prepare(db, select_sql.c_str(), -1, &stmt, NULL);
      sqlite3_step(stmt);

      while (sqlite3_column_type( stmt, 0 ) != SQLITE_NULL)
      {
        indices.push_back(sqlite3_column_int(stmt, 0));
        sqlite3_step(stmt);
      }

      // free memory and close connection
      sqlite3_finalize(stmt);
      sqlite3_close(db);

      return indices;
    }

    sqlite3* MzMLSqliteSwathHandler::openDB()
    {
      sqlite3 *db;
      int rc;

      // Open database
      rc = sqlite3_open(filename_.c_str(), &db);
      if (rc)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Can't open database: ") + sqlite3_errmsg(db));
      }
      return db;
    }

  } // namespace Internal
} // namespace OpenMS

