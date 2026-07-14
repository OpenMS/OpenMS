// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>
#include <vector>

namespace OpenMS
{
  /**
    @brief File adapter for Sqlite files

    This class contains certain helper functions to deal with Sqlite files.

    @note This is the public, SQLite-free interface. The raw SQLite C API
    (`sqlite3*`/`sqlite3_stmt*` handling and the @c Internal::SqliteHelper
    statement helpers) lives in the non-installed implementation header
    <OpenMS/FORMAT/SqliteConnector_impl.h>, which must only be included from
    .cpp files inside libOpenMS. Keeping the SQLite types out of this header is
    what allows SQLite to be a fully private dependency of OpenMS.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI SqliteConnector
  {
  public:

    /// how an sqlite db should be opened
    enum class SqlOpenMode
    {
      READ_ONLY,  ///< the DB must exist and is read-only
#ifndef READONLY
      READONLY = READ_ONLY, ///< compatibility alias for READ_ONLY
#endif
      READWRITE, ///< the DB is readable and writable, but must exist when opening it
      READWRITE_OR_CREATE ///< the DB readable and writable and is created new if not present already
    };

    /// Default constructor
    SqliteConnector() = delete;

    /// Constructor which opens a connection to @p filename
    /// @throws Exception::SqlOperationFailed if the file does not exist/cannot be created (depending on @p mode)
    explicit SqliteConnector(const std::string& filename, const SqlOpenMode mode = SqlOpenMode::READWRITE_OR_CREATE);

    /// Destructor
    ~SqliteConnector();

    /// @cond INTERNAL
    /**
      @brief Returns the raw native SQLite handle as an opaque pointer.

      The returned pointer is really an @c sqlite3*; it is exposed here as a
      @c void* so that this installed header stays free of any SQLite type.
      To obtain the typed handle, include the non-installed
      <OpenMS/FORMAT/SqliteConnector_impl.h> and call
      @c Internal::SqliteHelper::getNativeHandle(*this).

      @note The handle is tied to the lifetime of the SqliteConnector object;
      do not use it after the object has gone out of scope!
    */
    void* nativeHandle() const
    {
      return db_;
    }
    /// @endcond

    /**
      @brief Checks whether the given table exists

      @p tablename The name of the table to be checked

      @returns Whether the table exists or not
    */
    bool tableExists(const std::string& tablename);

    /// Counts the number of entries in SQL table @p table_name
    /// @throws Exception::SqlOperationFailed if table is unknown
    Size countTableRows(const std::string& table_name);

    /**
      @brief Checks whether the given table contains a certain column

      @p tablename The name of the table (needs to exist)
      @p colname The name of the column to be checked

      @returns Whether the column exists or not
    */
    bool columnExists(const std::string& tablename, const std::string& colname);

    /**
      @brief Executes a given SQL statement (insert statement)

      This is useful for writing a single row of data

      @p statement The SQL statement

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    void executeStatement(const std::string& statement);

    /**
      @brief Executes raw data SQL statements (insert statements)

      This is useful for a case where raw data should be inserted into sqlite
      databases, and the raw data needs to be passed separately as it cannot be
      part of a true SQL statement

        INSERT INTO TBL (ID, DATA) VALUES (100, ?1), (101, ?2), (102, ?3)"

      See also https://www.sqlite.org/c3ref/bind_blob.html

      @p prepare_statement The SQL statement
      @p data The data to bind

      @exception Exception::IllegalArgument is thrown if the SQL command fails.
    */
    void executeBindStatement(const std::string& prepare_statement, const std::vector<std::string>& data);

  protected:

    /**
      @brief Opens a new SQLite database

      @param[in] filename Filename of the database
      @param[in] mode See SqlOpenMode

      @note Call this only once!
    */
    void openDatabase_(const std::string& filename, const SqlOpenMode mode);

    void* db_ = nullptr; ///< opaque native SQLite handle (really an sqlite3*); see nativeHandle()

  };

} // namespace OpenMS
