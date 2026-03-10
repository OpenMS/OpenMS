// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/SYSTEM/File.h>

#include <sqlite3.h>

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathOSWWriter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OpenSwathOSWWriter* ptr = nullptr;
OpenSwathOSWWriter* nullPointer = nullptr;

START_SECTION(OpenSwathOSWWriter(const String& output_filename, bool uis_scores))
{
  ptr = new OpenSwathOSWWriter("", false);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~OpenSwathOSWWriter())
{
  delete ptr;
}
END_SECTION

START_SECTION(bool isActive() const)
{
  OpenSwathOSWWriter inactive_writer("", false);
  TEST_EQUAL(inactive_writer.isActive(), false)

  String temp_file = File::getTemporaryFile();
  OpenSwathOSWWriter active_writer(temp_file, false);
  TEST_EQUAL(active_writer.isActive(), true)
  File::remove(temp_file);
}
END_SECTION

START_SECTION([EXTRA] RUN.ID is stored as INTEGER not BLOB (regression test for BLOB/INTEGER join bug))
{
  // Regression test: RUN.ID was previously stored as BLOB because addRun()
  // passed the ID as a String to executeBindStatement(), which used
  // sqlite3_bind_blob for all parameters. FEATURE.RUN_ID was correctly stored
  // as INTEGER (embedded as a literal in SQL text by prepareLine()). This type
  // mismatch caused JOIN queries on FEATURE.RUN_ID = RUN.ID to return 0 rows.
  //
  // The fix: embed the integer ID as a literal in the SQL text in addRun(),
  // consistent with how prepareLine() writes FEATURE.RUN_ID.
  //
  // This test verifies:
  //   1. typeof(RUN.ID) == 'integer'     (was 'blob' before the fix)
  //   2. A JOIN between RUN and FEATURE on run_id returns rows when expected

  String temp_file = File::getTemporaryFile();

  // Use a large 64-bit value that previously triggered the BLOB storage bug
  const UInt64 large_run_id = 6130996817540441879ULL;

  // --- Write ---
  {
    OpenSwathOSWWriter writer(temp_file, false);
    writer.writeHeader();
    writer.addRun(large_run_id, "test_input_file.mzML");

    // Manually insert a minimal FEATURE row referencing the same run_id.
    // clearSignBit masks the sign bit so the stored id may be slightly smaller.
    // We use the same masking logic here.
    const UInt64 rid = large_run_id & ~(1ULL << 63); // mirrors clearSignBit
    SqliteConnector conn(temp_file);
    conn.executeStatement(
      "INSERT INTO FEATURE (ID, RUN_ID, PRECURSOR_ID, EXP_RT, EXP_IM, NORM_RT, DELTA_RT, LEFT_WIDTH, RIGHT_WIDTH, EXP_IM_LEFTWIDTH, EXP_IM_RIGHTWIDTH) "
      "VALUES (" + String(rid + 1) + ", " + String(rid) + ", 999, 100.0, NULL, 100.0, 0.0, 90.0, 110.0, NULL, NULL);"
    );
  }

  // --- Verify ---
  {
    SqliteConnector conn(temp_file, SqliteConnector::SqlOpenMode::READONLY);
    sqlite3* db = conn.getDB();

    // 1. Check typeof(RUN.ID) is 'integer'
    {
      sqlite3_stmt* stmt = nullptr;
      SqliteConnector::prepareStatement(db, &stmt, "SELECT typeof(ID) FROM RUN LIMIT 1;");
      sqlite3_step(stmt);
      const char* type_str = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
      TEST_NOT_EQUAL(type_str, nullptr)
      String type_value(type_str);
      sqlite3_finalize(stmt);
      TEST_EQUAL(type_value, "integer")
    }

    // 2. Check that JOIN on FEATURE.RUN_ID = RUN.ID returns rows
    {
      sqlite3_stmt* stmt = nullptr;
      SqliteConnector::prepareStatement(db, &stmt,
        "SELECT COUNT(*) FROM RUN INNER JOIN FEATURE ON FEATURE.RUN_ID = RUN.ID;");
      sqlite3_step(stmt);
      int join_count = sqlite3_column_int(stmt, 0);
      sqlite3_finalize(stmt);
      TEST_EQUAL(join_count, 1)
    }

    // 3. Verify the stored RUN.ID value round-trips correctly
    {
      const UInt64 rid = large_run_id & ~(1ULL << 63);
      sqlite3_stmt* stmt = nullptr;
      SqliteConnector::prepareStatement(db, &stmt, "SELECT ID FROM RUN LIMIT 1;");
      sqlite3_step(stmt);
      int64_t stored_id = sqlite3_column_int64(stmt, 0);
      sqlite3_finalize(stmt);
      TEST_EQUAL(static_cast<UInt64>(stored_id), rid)
    }
  }

  File::remove(temp_file);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
