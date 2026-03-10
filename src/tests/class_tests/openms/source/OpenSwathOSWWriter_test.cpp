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
  // Use READWRITE to allow creating temp tables for verification queries
  {
    SqliteConnector conn(temp_file);

    // 1. Check typeof(RUN.ID) is 'integer' (was 'blob' before the fix)
    conn.executeStatement(
      "CREATE TEMP TABLE _typeof_ok AS SELECT 1 FROM RUN WHERE typeof(ID) = 'integer';");
    TEST_EQUAL(conn.countTableRows("_typeof_ok"), 1)

    // 2. Check that JOIN on FEATURE.RUN_ID = RUN.ID returns rows
    conn.executeStatement(
      "CREATE TEMP TABLE _join_result AS "
      "SELECT RUN.ID FROM RUN INNER JOIN FEATURE ON FEATURE.RUN_ID = RUN.ID;");
    TEST_EQUAL(conn.countTableRows("_join_result"), 1)

    // 3. Verify the stored RUN.ID value round-trips correctly
    const UInt64 rid = large_run_id & ~(1ULL << 63);
    conn.executeStatement(
      "CREATE TEMP TABLE _id_match AS SELECT 1 FROM RUN WHERE ID = " + String(rid) + ";");
    TEST_EQUAL(conn.countTableRows("_id_match"), 1)
  }

  File::remove(temp_file);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
