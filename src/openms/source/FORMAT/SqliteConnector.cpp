// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/SqliteConnector.h>

#include <OpenMS/CONCEPT/Exception.h>

#include <sqlite3.h>

#include <cstring> // for strcmp
#include <iostream>

namespace OpenMS
{
  SqliteConnector::SqliteConnector(const String& filename, const SqlOpenMode mode)
  {
    openDatabase_(filename, mode);
  }
  SqliteConnector::~SqliteConnector()
  {
    int rc = sqlite3_close_v2(db_);
    if (rc != SQLITE_OK)
    {
      std::cout << " Encountered error in ~SqliteConnector: " << rc << std::endl;
    }
  }

  void SqliteConnector::openDatabase_(const String& filename, const SqlOpenMode mode)
  {
    // Open database
    int flags = 0;
    switch (mode)
    {
      case SqlOpenMode::READONLY:
        flags = SQLITE_OPEN_READONLY;
        break;
      case SqlOpenMode::READWRITE:
        flags = SQLITE_OPEN_READWRITE;
        break;
      case SqlOpenMode::READWRITE_OR_CREATE:
        flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE;
        break;
    }
    int rc = sqlite3_open_v2(filename.c_str(), &db_, flags, nullptr);
    if (rc)
    {
      throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not open sqlite db '" + filename + "' in mode " + String(int(mode)));
    }
  }

  bool SqliteConnector::columnExists(sqlite3 *db, const String& tablename, const String& colname)
  {
    bool found = false;

    sqlite3_stmt* xcntstmt;
    prepareStatement(db, &xcntstmt, "PRAGMA table_info(" + tablename + ")");

    // Go through all columns and check whether the required column exists
    sqlite3_step(xcntstmt);
    while (sqlite3_column_type(xcntstmt, 0) != SQLITE_NULL)
    {
      if (strcmp(colname.c_str(), reinterpret_cast<const char*>(sqlite3_column_text(xcntstmt, 1))) == 0)
      {
        found = true;
        break;
      }
      sqlite3_step(xcntstmt);
    }
    sqlite3_finalize(xcntstmt);

    return found;
  }

  bool SqliteConnector::tableExists(sqlite3 *db, const String& tablename)
  {
    sqlite3_stmt* stmt;
    prepareStatement(db, &stmt, "SELECT 1 FROM sqlite_master WHERE type='table' AND name='" + tablename + "';");

    sqlite3_step(stmt);
    // if we get a row back, the table exists:
    bool found = (sqlite3_column_type(stmt, 0) != SQLITE_NULL);
    sqlite3_finalize(stmt);

    return found;
  }

  Size SqliteConnector::countTableRows(const String& table_name)
  {
    sqlite3_stmt* stmt;
    String select_runs = "SELECT count(*) FROM " + table_name + ";";
    this->prepareStatement(&stmt, select_runs);
    sqlite3_step(stmt);
    if (sqlite3_column_type(stmt, 0) == SQLITE_NULL)
    {
      throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Could not retrieve " + table_name + " table count!");
    }
    Size res = sqlite3_column_int64(stmt, 0);
    sqlite3_finalize(stmt);
    return res;
  }

  void SqliteConnector::executeStatement(sqlite3 *db, const String& statement)
  {
    char *zErrMsg = nullptr;
    int rc = sqlite3_exec(db, statement.c_str(), nullptr /* callback */, nullptr, &zErrMsg);
    if (rc != SQLITE_OK)
    {
      String error(zErrMsg);
      std::cerr << "Error message after sqlite3_exec" << std::endl;
      std::cerr << "Prepared statement " << statement << std::endl;
      sqlite3_free(zErrMsg);
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
    }
  }

  void SqliteConnector::prepareStatement(sqlite3 *db, sqlite3_stmt** stmt, const String& prepare_statement)
  {
    int rc = sqlite3_prepare_v2(db, prepare_statement.c_str(), (int)prepare_statement.size(), stmt, nullptr);
    if (rc != SQLITE_OK)
    {
      std::cerr << "Error message after sqlite3_prepare_v2" << std::endl;
      std::cerr << "Prepared statement " << prepare_statement << std::endl;
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
    }
  }

  void SqliteConnector::executeBindStatement(sqlite3 *db, const String& prepare_statement, const std::vector<String>& data)
  {
    int rc;
    sqlite3_stmt *stmt = nullptr;
    prepareStatement(db, &stmt, prepare_statement);
    for (Size k = 0; k < data.size(); k++)
    {
      // Fifth argument is a destructor for the blob.
      // SQLITE_STATIC because the statement is finalized
      // before the buffer is freed:
      rc = sqlite3_bind_blob(stmt, k+1, data[k].c_str(), (int)data[k].size(), SQLITE_STATIC);
      if (rc != SQLITE_OK)
      {
        std::cerr << "SQL error after sqlite3_bind_blob at iteration " << k << std::endl;
        std::cerr << "Prepared statement " << prepare_statement << std::endl;
        // TODO this is a mem-leak (missing sqlite3_finalize())
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
      }
    }

    rc = sqlite3_step(stmt);
    if (rc != SQLITE_DONE)
    {
      std::cerr << "SQL error after sqlite3_step" << std::endl;
      std::cerr << "Prepared statement " << prepare_statement << std::endl;
      // TODO this is a mem-leak (missing sqlite3_finalize())
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, sqlite3_errmsg(db));
    }

    // free memory
    sqlite3_finalize(stmt);
  }

  namespace Internal::SqliteHelper
    {

      template <> bool extractValue<double>(double* dst, sqlite3_stmt* stmt, int pos) //explicit specialization
      {
        if (sqlite3_column_type(stmt, pos) != SQLITE_NULL)
        {
          *dst = sqlite3_column_double(stmt, pos);
          return true;
        }
        return false;
      }

      template <> bool extractValue<int>(Int32* dst, sqlite3_stmt* stmt, int pos) //explicit specialization
      {
        if (sqlite3_column_type(stmt, pos) != SQLITE_NULL)
        {
          *dst = sqlite3_column_int(stmt, pos); // sqlite3_column_int returns 32bit integers
          return true;
        }
        return false;
      }
      template <> bool extractValue<Int64>(Int64* dst, sqlite3_stmt* stmt, int pos) //explicit specialization
      {
        if (sqlite3_column_type(stmt, pos) != SQLITE_NULL)
        {
          *dst = sqlite3_column_int64(stmt, pos);
          return true;
        }
        return false;
      }

      template <> bool extractValue<String>(String* dst, sqlite3_stmt* stmt, int pos) //explicit specialization
      {
        if (sqlite3_column_type(stmt, pos) != SQLITE_NULL)
        {
          *dst = String(reinterpret_cast<const char*>(sqlite3_column_text(stmt, pos)));
          return true;
        }
        return false;
      }

      template <> bool extractValue<std::string>(std::string* dst, sqlite3_stmt* stmt, int pos) //explicit specialization
      {
        if (sqlite3_column_type(stmt, pos) != SQLITE_NULL)
        {
          *dst = std::string(reinterpret_cast<const char*>(sqlite3_column_text(stmt, pos)));
          return true;
        }

        return false;
      }

      SqlState nextRow(sqlite3_stmt* stmt, SqlState current)
      {
        if (current != SqlState::SQL_ROW)
        { // querying a new row after the last invocation gave 'SQL_DONE' might loop around
          // to the first entry and give an infinite loop!!!
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sql operation requested on SQL_DONE/SQL_ERROR state. This should never happen. Please file a bug report!");
        }
        int rc = sqlite3_step(stmt);
        if (rc == SQLITE_ROW)
        {
          return SqlState::SQL_ROW;
        }
        if (rc == SQLITE_DONE)
        {
          return SqlState::SQL_DONE;
        }
        if (rc == SQLITE_ERROR)
        {
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sql operation failed with SQLITE_ERROR!");
        }
        if (rc == SQLITE_BUSY)
        {
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sql operation failed with SQLITE_BUSY!");
        }
        if (rc == SQLITE_MISUSE)
        {
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sql operation failed with SQLITE_MISUSE!");
        }
        throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Sql operation failed with unexpected error code!");
      }

      /// Special case: store integer in a string data value
      bool extractValueIntStr(String* dst, sqlite3_stmt* stmt, int pos)
      {
        if (sqlite3_column_type(stmt, pos) == SQLITE_INTEGER)
        {
          *dst = sqlite3_column_int(stmt, pos);
          return true;
        }
        return false;
      }

      double extractDouble(sqlite3_stmt* stmt, int pos)
      {
        double res;
        if (!extractValue<double>(&res, stmt, pos)) 
        {
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Conversion of column " + String(pos) + " to double failed");
        }
        return res;
      }

      float extractFloat(sqlite3_stmt* stmt, int pos)
      {
        double res; // there is no sqlite3_column_float.. so we extract double and convert
        if (!extractValue<double>(&res, stmt, pos))
        {
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Conversion of column " + String(pos) + " to double/float failed");
        }
        return (float)res;
      }

      int extractInt(sqlite3_stmt* stmt, int pos)
      {
        int res;
        if (!extractValue<int>(&res, stmt, pos))
        {
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Conversion of column " + String(pos) + " to int failed");
        }
        return res;
      }

      Int64 extractInt64(sqlite3_stmt* stmt, int pos)
      {
        Int64 res;
        if (!extractValue<Int64>(&res, stmt, pos))
        {
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Conversion of column " + String(pos) + " to Int64 failed");
        }
        return res;
      }

      String extractString(sqlite3_stmt* stmt, int pos)
      {
        String res;
        if (!extractValue<String>(&res, stmt, pos))
        {
          throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Conversion of column " + String(pos) + " to String failed");
        }
        return res;
      }

      char extractChar(sqlite3_stmt* stmt, int pos)
      {
        return extractString(stmt, pos)[0];
      }

      bool extractBool(sqlite3_stmt* stmt, int pos)
      {
        return extractInt(stmt, pos) != 0;
      }

    }

}


