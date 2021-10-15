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

#include <sqlite3.h>

namespace OpenMS
{

  namespace Internal
  {
    namespace SqliteHelper
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

      template <> bool extractValue<int>(int* dst, sqlite3_stmt* stmt, int pos) //explicit specialization
      {
        if (sqlite3_column_type(stmt, pos) != SQLITE_NULL)
        {
          *dst = sqlite3_column_int(stmt, pos);
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
} // namespace OpenMS
