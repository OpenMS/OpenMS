// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

// -----------------------------------------------------------------------------
// PRIVATE (non-installed) implementation header.
//
// This header exposes the raw SQLite C API (sqlite3 / sqlite3_stmt) that backs
// OpenMS::SqliteConnector. It is deliberately NOT installed, NOT exported and
// NOT part of the public API, so that SQLite can remain a fully private link
// dependency of libOpenMS (no SQLite type ever appears in an installed header).
//
// It must only be included from .cpp files that live inside libOpenMS. Public
// (installed) headers must never include it.
// -----------------------------------------------------------------------------

#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <sqlite3.h>

#include <sstream>
#include <string>
#include <type_traits> // for is_same
#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    namespace SqliteHelper
    {
      /// Retrieve the typed native SQLite handle behind a SqliteConnector.
      inline sqlite3* getNativeHandle(SqliteConnector& conn)
      {
        return static_cast<sqlite3*>(conn.nativeHandle());
      }

      /**
        @brief Checks whether the given table exists

        @p db The sqlite database (needs to be open)
        @p tablename The name of the table to be checked

        @returns Whether the table exists or not
      */
      bool tableExists(sqlite3* db, const std::string& tablename);

      /**
        @brief Checks whether the given table contains a certain column

        @p db The sqlite database (needs to be open)
        @p tablename The name of the table (needs to exist)
        @p colname The name of the column to be checked

        @returns Whether the column exists or not
      */
      bool columnExists(sqlite3* db, const std::string& tablename, const std::string& colname);

      /**
        @brief Executes a given SQL statement (insert statement)

        This wraps sqlite3_exec with proper error handling.

        @p db The sqlite database (needs to be open)
        @p statement The SQL statement

        @exception Exception::IllegalArgument is thrown if the SQL command fails.
      */
      void executeStatement(sqlite3* db, const std::string& statement);

      /**
        @brief Converts an SQL statement into a prepared statement

        Internally calls sqlite3_prepare_v2.

        @p db The sqlite database (needs to be open)
        @p stmt The prepared statement (output)
        @p prepare_statement The SQL statement to prepare (input)

        @exception Exception::IllegalArgument is thrown if the SQL command fails.
      */
      void prepareStatement(sqlite3* db, sqlite3_stmt** stmt, const std::string& prepare_statement);

      /// Convenience overload: prepare a statement directly on a SqliteConnector's native handle.
      inline void prepareStatement(SqliteConnector& conn, sqlite3_stmt** stmt, const std::string& prepare_statement)
      {
        prepareStatement(getNativeHandle(conn), stmt, prepare_statement);
      }

      /**
        @brief Executes raw data SQL statements (insert statements)

        See also https://www.sqlite.org/c3ref/bind_blob.html

        @p db The sqlite database (needs to be open)
        @p prepare_statement The SQL statement
        @p data The data to bind

        @exception Exception::IllegalArgument is thrown if the SQL command fails.
      */
      void executeBindStatement(sqlite3* db, const std::string& prepare_statement, const std::vector<std::string>& data);

      /// Sql only stores signed 64bit ints, so we remove the highest bit, because some/most
      /// of our sql-insert routines first convert to string, which might yield an uint64 which cannot
      /// be represented as int64, and sqlite would attempt to store it as double(!), which will loose precision
      template <typename T>
      UInt64 clearSignBit(T /*value*/)
      {
        static_assert(std::is_same<T, std::false_type>::value, "Wrong input type to clearSignBit(). Please pass unsigned 64bit ints!");
        return 0;
      };
      /// only allow UInt64 specialization
      template <>
      inline UInt64 clearSignBit(UInt64 value) {
        return value & ~(1ULL << 63);
      }


      enum class SqlState
      {
        SQL_ROW,
        SQL_DONE,
        SQL_ERROR ///< includes SQLITE_BUSY, SQLITE_ERROR, SQLITE_MISUSE
      };

      /**
        @brief retrieves the next row from a prepared statement

        If you receive 'SqlState::SQL_DONE', do NOT query nextRow() again,
        because you might enter an infinite loop!
        To avoid oversights, you can pass the old return value into the function again
        and get an Exception which will tell you that there is buggy code!

        @param[in] stmt Sqlite statement object
        @param[in] current Return value of the previous call to this function.
        @return one of SqlState::SQL_ROW or SqlState::SQL_DONE
        @throws Exception::SqlOperationFailed if state would be SqlState::ERROR
      */
      SqlState nextRow(sqlite3_stmt* stmt, SqlState current = SqlState::SQL_ROW);


      /**
        @brief Extracts a specific value from an SQL column

        @p dst Destination (where to store the value)
        @p stmt Sqlite statement object
        @p pos Column position

        For example, to extract a specific integer from column 5 of an SQL statement, one can use:

          sqlite3_stmt* stmt;
          sqlite3* db;
          SqliteHelper::prepareStatement(db, &stmt, select_sql);
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
      template <> bool extractValue<Int64>(Int64* dst, sqlite3_stmt* stmt, int pos); //explicit specialization

      template <> bool extractValue<std::string>(std::string* dst, sqlite3_stmt* stmt, int pos); //explicit specialization

      /// Special case where an integer should be stored in a std::string field
      bool extractValueIntStr(std::string* dst, sqlite3_stmt* stmt, int pos);

      /** @defgroup sqlThrowingGetters Functions for getting values from sql-select statements

          All these function throw Exception::SqlOperationFailed if the given position is of the wrong type.
       @{
       */
      double extractDouble(sqlite3_stmt* stmt, int pos);
      float extractFloat(sqlite3_stmt* stmt, int pos); ///< convenience function; note: in SQL there is no float, just double. So this might be narrowing.
      int extractInt(sqlite3_stmt* stmt, int pos);
      Int64 extractInt64(sqlite3_stmt* stmt, int pos);
      std::string extractString(sqlite3_stmt* stmt, int pos);
      char extractChar(sqlite3_stmt* stmt, int pos);
      bool extractBool(sqlite3_stmt* stmt, int pos);
      /** @} */ // end of sqlThrowingGetters
    }
  }
} // namespace OpenMS
