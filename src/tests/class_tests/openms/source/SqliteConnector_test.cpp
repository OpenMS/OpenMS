// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

// NOTE: this test deliberately includes only the public, SQLite-free header.
// SQLite is a private dependency of OpenMS, so downstream consumers (like this
// test) do not have <sqlite3.h> on their include path. The raw SQLite C API
// (Internal::SqliteHelper, sqlite3*/sqlite3_stmt*) lives in the non-installed
// <OpenMS/FORMAT/SqliteConnector_impl.h> and is exercised transitively by the
// OSWFile / MzMLSqlite / TransitionPQP tests. Here we cover the public API only.
#include <OpenMS/FORMAT/SqliteConnector.h>

#include <vector>
#include <string>

///////////////////////////

START_TEST(SqliteConnector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

using Mode = SqliteConnector::SqlOpenMode;

// temporary database filename (created lazily by SqliteConnector in CREATE mode)
std::string tmp_db;
NEW_TMP_FILE(tmp_db)

SqliteConnector* ptr = nullptr;
SqliteConnector* null_ptr = nullptr;

START_SECTION((SqliteConnector(const std::string& filename, const SqlOpenMode mode = SqlOpenMode::READWRITE_OR_CREATE)))
{
  // READWRITE_OR_CREATE creates the database if it does not exist
  ptr = new SqliteConnector(tmp_db, Mode::READWRITE_OR_CREATE);
  TEST_NOT_EQUAL(ptr, null_ptr)
  // a successfully opened connection is usable (the native handle is fully private)
  ptr->executeStatement("CREATE TABLE IF NOT EXISTS liveness_check (ID INT)");
  TEST_EQUAL(ptr->tableExists("liveness_check"), true)
}
END_SECTION

START_SECTION((~SqliteConnector()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void executeStatement(const std::string& statement)))
{
  SqliteConnector conn(tmp_db, Mode::READWRITE_OR_CREATE);
  conn.executeStatement("DROP TABLE IF EXISTS T");
  conn.executeStatement("CREATE TABLE T (ID INT, NAME TEXT)");
  conn.executeStatement("INSERT INTO T (ID, NAME) VALUES (1, 'a')");
  TEST_EQUAL(conn.tableExists("T"), true)

  // malformed SQL -> IllegalArgument (exercises the internal prepare/exec error path)
  TEST_EXCEPTION(Exception::IllegalArgument, conn.executeStatement("THIS IS NOT SQL"))
}
END_SECTION

START_SECTION((bool tableExists(const std::string& tablename)))
{
  SqliteConnector conn(tmp_db, Mode::READWRITE_OR_CREATE);
  conn.executeStatement("DROP TABLE IF EXISTS T");
  conn.executeStatement("CREATE TABLE T (ID INT, NAME TEXT)");
  TEST_EQUAL(conn.tableExists("T"), true)
  TEST_EQUAL(conn.tableExists("DOES_NOT_EXIST"), false)
}
END_SECTION

START_SECTION((bool columnExists(const std::string& tablename, const std::string& colname)))
{
  SqliteConnector conn(tmp_db, Mode::READWRITE_OR_CREATE);
  conn.executeStatement("DROP TABLE IF EXISTS T");
  conn.executeStatement("CREATE TABLE T (ID INT, NAME TEXT)");
  TEST_EQUAL(conn.columnExists("T", "ID"), true)
  TEST_EQUAL(conn.columnExists("T", "NAME"), true)
  TEST_EQUAL(conn.columnExists("T", "MISSING"), false)
}
END_SECTION

START_SECTION((Size countTableRows(const std::string& table_name)))
{
  SqliteConnector conn(tmp_db, Mode::READWRITE_OR_CREATE);
  conn.executeStatement("DROP TABLE IF EXISTS T");
  conn.executeStatement("CREATE TABLE T (ID INT, NAME TEXT)");
  TEST_EQUAL(conn.countTableRows("T"), 0)
  conn.executeStatement("INSERT INTO T (ID, NAME) VALUES (1, 'a')");
  conn.executeStatement("INSERT INTO T (ID, NAME) VALUES (2, 'b')");
  TEST_EQUAL(conn.countTableRows("T"), 2)

  // unknown table cannot be prepared -> IllegalArgument (exercises the internal prepare error path)
  TEST_EXCEPTION(Exception::IllegalArgument, conn.countTableRows("UNKNOWN"))
}
END_SECTION

START_SECTION((void executeBindStatement(const std::string& prepare_statement, const std::vector<std::string>& data)))
{
  SqliteConnector conn(tmp_db, Mode::READWRITE_OR_CREATE);
  conn.executeStatement("DROP TABLE IF EXISTS T");
  conn.executeStatement("CREATE TABLE T (ID INT, NAME TEXT)");
  conn.executeBindStatement("INSERT INTO T (ID, NAME) VALUES (1, ?1)", std::vector<std::string>{"bound_value"});
  TEST_EQUAL(conn.countTableRows("T"), 1)
}
END_SECTION

START_SECTION([EXTRA] (open modes: existing DB read-only and missing-file errors))
{
  // Populate a database, then reopen it read-only and read back
  {
    SqliteConnector conn(tmp_db, Mode::READWRITE_OR_CREATE);
    conn.executeStatement("DROP TABLE IF EXISTS RO");
    conn.executeStatement("CREATE TABLE RO (ID INT)");
    conn.executeStatement("INSERT INTO RO (ID) VALUES (1)");
    conn.executeStatement("INSERT INTO RO (ID) VALUES (2)");
    conn.executeStatement("INSERT INTO RO (ID) VALUES (3)");
  }
  {
    SqliteConnector conn(tmp_db, Mode::READ_ONLY);
    TEST_EQUAL(conn.tableExists("RO"), true)
    TEST_EQUAL(conn.countTableRows("RO"), 3)
  }

  // Opening a non-existent database without CREATE must fail
  std::string missing;
  NEW_TMP_FILE(missing) // a unique filename that is never created
  TEST_EXCEPTION(Exception::SqlOperationFailed, SqliteConnector(missing, Mode::READ_ONLY))
  TEST_EXCEPTION(Exception::SqlOperationFailed, SqliteConnector(missing, Mode::READWRITE))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
