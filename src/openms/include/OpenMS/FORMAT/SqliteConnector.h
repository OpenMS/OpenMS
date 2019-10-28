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
    ~SqliteConnector();

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
      @brief Checks whether the given table exists

      @p tablename The name of the table to be checked

      @returns Whether the table exists or not
    */
    bool tableExists(const String& tablename)
    {
      return tableExists(db_, tablename);
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
    void prepareStatement(sqlite3_stmt** stmt, const String& prepare_statement)
    {
      prepareStatement(db_, stmt, prepare_statement);
    }

    /**
      @brief Checks whether the given table exists

      @p db The sqlite database (needs to be open)
      @p tablename The name of the table to be checked

      @returns Whether the table exists or not
    */
    static bool tableExists(sqlite3* db, const String& tablename);

    /**
      @brief Checks whether the given table contains a certain column

      @p db The sqlite database (needs to be open)
      @p tablename The name of the table (needs to exist)
      @p colname The name of the column to be checked

      @returns Whether the column exists or not
    */
    static bool columnExists(sqlite3* db, const String& tablename, const String& colname);

    /**
      @brief Executes a given SQL statement (insert statement)

      This is useful for writing a single row of data. It wraps sqlite3_exec with proper error handling.

      @p db The sqlite database (needs to be open)
      @p statement The SQL statement

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    static void executeStatement(sqlite3* db, const std::stringstream& statement);

    /**
      @brief Executes a given SQL statement (insert statement)

      This is useful for writing a single row of data. It wraps sqlite3_exec with proper error handling.

      @p db The sqlite database (needs to be open)
      @p statement The SQL statement

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    static void executeStatement(sqlite3* db, const String& statement);

    /**
      @brief Converts an SQL statement into a prepared statement

      This routine converts SQL text into a prepared statement object and
      returns a pointer to that object. This interface requires a database
      connection created by a prior call to sqlite3_open() and a text string
      containing the SQL statement to be prepared. This API does not actually
      evaluate the SQL statement. It merely prepares the SQL statement for
      evaluation.

      This is useful for handling errors in a consistent manner. Internally
      calls sqlite3_prepare_v2.

      @p db The sqlite database (needs to be open)
      @p stmt The prepared statement (output)
      @p prepare_statement The SQL statement to prepare (input)

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    static void prepareStatement(sqlite3* db, sqlite3_stmt** stmt, const String& prepare_statement);


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
    static void executeBindStatement(sqlite3* db, const String& prepare_statement, const std::vector<String>& data);

protected:

    /**
      @brief Opens a new SQLite database

      @filename Filename of the database

      @note Call this only once!
    */
    void openDatabase(const String& filename);

protected:
    sqlite3 *db_;

  };

  namespace Internal
  {
    namespace SqliteHelper
    {

      /**
        @brief Extracts a specific value from an SQL column

        @p dst Destination (where to store the value)
        @p stmt Sqlite statement object
        @p pos Column position

        For example, to extract a specific integer from column 5 of an SQL statement, one can use:

          sqlite3_stmt* stmt;
          sqlite3* db;
          SqliteConnector::prepareStatement(db, &stmt, select_sql);
          sqlite3_step(stmt);

          double target;
          while (sqlite3_column_type(stmt, 0) != SQLITE_NULL)
          {
            extractValue<double>(&target, stmt, 5);
            sqlite3_step( stmt );
          }
          sqlite3_finalize(stmt);
      */
      template <typename ValueType>
      bool extractValue(ValueType* /* dst */, sqlite3_stmt* /* stmt */, int /* pos */)
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Not implemented");
      }

      template <> bool extractValue<double>(double* dst, sqlite3_stmt* stmt, int pos); //explicit specialization

      template <> bool extractValue<int>(int* dst, sqlite3_stmt* stmt, int pos); //explicit specialization

      template <> bool extractValue<String>(String* dst, sqlite3_stmt* stmt, int pos); //explicit specialization

      template <> bool extractValue<std::string>(std::string* dst, sqlite3_stmt* stmt, int pos); //explicit specialization

      /// Special case where an integer should be stored in a String field
      bool extractValueIntStr(String* dst, sqlite3_stmt* stmt, int pos);

    }
  }


} // namespace OpenMS


