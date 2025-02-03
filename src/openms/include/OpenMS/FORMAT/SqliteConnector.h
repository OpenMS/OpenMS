// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <type_traits> // for is_same

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

    /// how an sqlite db should be opened
    enum class SqlOpenMode
    {
      READONLY,  ///< the DB must exist and is read-only
      READWRITE, ///< the DB is readable and writable, but must exist when opening it
      READWRITE_OR_CREATE ///< the DB readable and writable and is created new if not present already
    };

    /// Default constructor
    SqliteConnector() = delete;

    /// Constructor which opens a connection to @p filename
    /// @throws Exception::SqlOperationFailed if the file does not exist/cannot be created (depending on @p mode)
    explicit SqliteConnector(const String& filename, const SqlOpenMode mode = SqlOpenMode::READWRITE_OR_CREATE);

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

    /// Counts the number of entries in SQL table @p table_name
    /// @throws Exception::SqlOperationFailed if table is unknown
    Size countTableRows(const String& table_name);

    /**
      @brief Checks whether the given table contains a certain column

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

      @param filename Filename of the database
      @param mode See SqlOpenMode

      @note Call this only once!
    */
    void openDatabase_(const String& filename, const SqlOpenMode mode);

  protected:
    sqlite3* db_ = nullptr;

  };

  namespace Internal
  {
    namespace SqliteHelper
    {
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

        @param stmt Sqlite statement object
        @param current Return value of the previous call to this function.
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
      template <> bool extractValue<Int64>(Int64* dst, sqlite3_stmt* stmt, int pos); //explicit specialization

      template <> bool extractValue<String>(String* dst, sqlite3_stmt* stmt, int pos); //explicit specialization

      template <> bool extractValue<std::string>(std::string* dst, sqlite3_stmt* stmt, int pos); //explicit specialization

      /// Special case where an integer should be stored in a String field
      bool extractValueIntStr(String* dst, sqlite3_stmt* stmt, int pos);

      /** @defgroup sqlThrowingGetters Functions for getting values from sql-select statements

          All these function throw Exception::SqlOperationFailed if the given position is of the wrong type.
       @{
       */
      double extractDouble(sqlite3_stmt* stmt, int pos);
      float extractFloat(sqlite3_stmt* stmt, int pos); ///< convenience function; note: in SQL there is no float, just double. So this might be narrowing.
      int extractInt(sqlite3_stmt* stmt, int pos);
      Int64 extractInt64(sqlite3_stmt* stmt, int pos);
      String extractString(sqlite3_stmt* stmt, int pos);
      char extractChar(sqlite3_stmt* stmt, int pos);
      bool extractBool(sqlite3_stmt* stmt, int pos);
      /** @} */ // end of sqlThrowingGetters
    }
  }


} // namespace OpenMS
