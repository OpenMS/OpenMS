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

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <iostream>
#include <sstream>

// forward declarations
struct sqlite3;
struct sqlite3_stmt;

#include <sqlite3.h>

namespace OpenMS
{
  /**
    @brief File adapter for Sqlite files

    This class contains certain helper functions to deal with Sqlite files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI SqliteConnector
  {
public:

    /// Default constructor
    SqliteConnector();

    explicit SqliteConnector(const String& filename)
    {
      openDatabase(filename);
    }

    /// Destructor
    ~SqliteConnector()
    {
      sqlite3_close(db_);
    }

    /**
      @brief Returns the raw pointer to the database 

      @note The pointer is tied to the lifetime of the SqliteConnector object,
      do not use it after the object has gone out of scope!

      @returns SQLite database ptr

    */
    sqlite3* getDB()
    {
      return db_;
    }

    /**
      @brief Checkes whether the given table contains a certain column

      @p tablename The name of the table (needs to exist)
      @p colname The name of the column to be checked

      @returns Whether the column exists or not
    */
    bool columnExists(const String& tablename, const String& colname)
    {
      return columnExists(db_, tablename, colname);
    }

    /**
      @brief Executes a given SQL statement (insert statement)

      This is useful for writing a single row of data

      @p statement The SQL statement

    */
    void executeStatement(const std::stringstream& statement)
    {
      executeStatement(db_, statement);
    }

    /**
      @brief Executes a given SQL statement (insert statement)

      This is useful for writing a single row of data

      @p statement The SQL statement

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    void executeStatement(const String& statement)
    {
      executeStatement(db_, statement);
    }

    /**
      @brief Executes raw data SQL statements (insert statements) 

      This is useful for a case where raw data should be inserted into sqlite
      databases, and the raw data needs to be passed separately as it cannot be
      part of a true SQL statement

        INSERT INTO TBL (ID, DATA) VALUES (100, ?1), (101, ?2), (102, ?3)" 

      See also https://www.sqlite.org/c3ref/bind_blob.html

      @p statement The SQL statement
      @p data The data to bind

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    void executeBindStatement(const String& prepare_statement, const std::vector<String>& data)
    {
      executeBindStatement(db_, prepare_statement, data);
    }

    /**
      @brief Prepares a SQL statement

      This is useful for handling errors in a consistent manner.

      @p db The sqlite database (needs to be open)
      @p statement The SQL statement
      @p data The data to bind

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    void executePreparedStatement(sqlite3_stmt *stmt, const String& prepare_statement)
    {
      executePreparedStatement(db_, stmt, prepare_statement);
    }




    /**
      @brief Checkes whether the given table contains a certain column

      @p db The sqlite database (needs to be open)
      @p tablename The name of the table (needs to exist)
      @p colname The name of the column to be checked

      @returns Whether the column exists or not
    */
    static bool columnExists(sqlite3 *db, const String& tablename, const String& colname);

    /**
      @brief Executes a given SQL statement (insert statement)

      This is useful for writing a single row of data

      @p db The sqlite database (needs to be open)
      @p statement The SQL statement

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    static void executeStatement(sqlite3 *db, const std::stringstream& statement)
    {
      char *zErrMsg = nullptr;
      std::string insert_str = statement.str();
      int rc = sqlite3_exec(db, insert_str.c_str(), nullptr /* callback */, nullptr, &zErrMsg);
      if (rc != SQLITE_OK)
      {
        String error (zErrMsg);
        std::cerr << "Error message after sqlite3_exec" << std::endl;
        std::cerr << "Prepared statement " << statement.str() << std::endl;
        sqlite3_free(zErrMsg);
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
      }
    }

    /**
      @brief Executes a given SQL statement (insert statement)

      This is useful for writing a single row of data. It wraps sqlite3_exec with proper error handling.

      @p db The sqlite database (needs to be open)
      @p statement The SQL statement

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    static void executeStatement(sqlite3 *db, const String& statement)
    {
      char *zErrMsg = nullptr;
      int rc = sqlite3_exec(db, statement.c_str(), nullptr /* callback */, nullptr, &zErrMsg);
      if (rc != SQLITE_OK)
      {
        String error (zErrMsg);
        std::cerr << "Error message after sqlite3_exec" << std::endl;
        std::cerr << "Prepared statement " << statement << std::endl;
        sqlite3_free(zErrMsg);
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
      }
    }

    /**
      @brief Prepares a SQL statement

      This is useful for handling errors in a consistent manner.

      @p db The sqlite database (needs to be open)
      @p statement The SQL statement
      @p data The data to bind

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    static void executePreparedStatement(sqlite3 *db, sqlite3_stmt *stmt, const String& prepare_statement)
    {
      int rc = sqlite3_prepare_v2(db, prepare_statement.c_str(), prepare_statement.size(), &stmt, nullptr);
      if (rc != SQLITE_OK)
      {
        std::cerr << "Error message after sqlite3_prepare_v2" << std::endl;
        std::cerr << "Prepared statement " << prepare_statement << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      }
    }

    /**
      @brief Executes raw data SQL statements (insert statements) 

      This is useful for a case where raw data should be inserted into sqlite
      databases, and the raw data needs to be passed separately as it cannot be
      part of a true SQL statement

        INSERT INTO TBL (ID, DATA) VALUES (100, ?1), (101, ?2), (102, ?3)" 

      See also https://www.sqlite.org/c3ref/bind_blob.html

      @p db The sqlite database (needs to be open)
      @p statement The SQL statement
      @p data The data to bind

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    static void executeBindStatement(sqlite3 *db, const String& prepare_statement, const std::vector<String>& data)
    {
      sqlite3_stmt *stmt = nullptr;
      const char *curr_loc;
      int rc = sqlite3_prepare_v2(db, prepare_statement.c_str(), prepare_statement.size(), &stmt, &curr_loc);
      if (rc != SQLITE_OK)
      {
        std::cerr << "Error message after sqlite3_prepare_v2" << std::endl;
        std::cerr << "Prepared statement " << prepare_statement << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      }

      for (Size k = 0; k < data.size(); k++)
      {
        // Fifth argument is a destructor for the blob.
        // SQLITE_STATIC because the statement is finalized
        // before the buffer is freed:
        rc = sqlite3_bind_blob(stmt, k+1, data[k].c_str(), data[k].size(), SQLITE_STATIC);
        if (rc != SQLITE_OK)
        {
          std::cerr << "SQL error after sqlite3_bind_blob at iteration " << k << std::endl;
          std::cerr << "Prepared statement " << prepare_statement << std::endl;
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
        } 
      }

      rc = sqlite3_step(stmt);
      if (rc != SQLITE_DONE)
      {
        std::cerr << "SQL error after sqlite3_step" << std::endl;
        std::cerr << "Prepared statement " << prepare_statement << std::endl;
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      }

      // free memory
      sqlite3_finalize(stmt);
    }

protected:

    /**
      @brief Opens a new SQLite database

      @filename Filename of the database

      @note Call this only once!
    */
    void openDatabase(const String& filename)
    {
      // Open database
      int rc = sqlite3_open(filename.c_str(), &db_);
      if (rc)
      {
        throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db_));
      }
    }

protected:
    sqlite3 *db_;

  };

} // namespace OpenMS


