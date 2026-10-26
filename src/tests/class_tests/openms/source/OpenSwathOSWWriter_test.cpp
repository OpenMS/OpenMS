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
#include <OpenMS/DATASTRUCTURES/DataValue.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/TempFiles.h>

#include <string>
#include <utility>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathOSWWriter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OpenSwathOSWWriter* ptr = nullptr;
OpenSwathOSWWriter* nullPointer = nullptr;

START_SECTION(OpenSwathOSWWriter(const std::string& output_filename, bool uis_scores))
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

  std::string temp_file = TempFiles::getTemporaryFile();
  OpenSwathOSWWriter active_writer(temp_file, false);
  TEST_EQUAL(active_writer.isActive(), true)
  File::remove(temp_file);
}
END_SECTION

START_SECTION([EXTRA] RUN.ID is stored as INTEGER not BLOB (regression test for BLOB/INTEGER join bug))
{
  // Regression test: RUN.ID was previously stored as BLOB because addRun()
  // passed the ID as a String to executeBindStatement(), which used
  // sqlite3_bind_blob for all parameters. That type mismatch caused JOIN
  // queries on FEATURE.RUN_ID = RUN.ID to return 0 rows.
  //
  // This test verifies:
  //   1. typeof(RUN.ID) == 'integer'     (was 'blob' before the fix)
  //   2. A JOIN between RUN and FEATURE on run_id returns rows when expected

  std::string temp_file = TempFiles::getTemporaryFile();

  // Use a large 64-bit value that previously triggered the BLOB storage bug
  const UInt64 large_run_id = 6130996817540441879ULL;

  // --- Write ---
  {
    OpenSwathOSWWriter writer(temp_file, false);
    writer.writeHeader();
    writer.addRun(large_run_id, "test_input_file.mzML");

    // Insert a minimal FEATURE row through the prepared-row writer path.
    // clearSignBit masks the sign bit so the stored id may be slightly smaller.
    // We use the same masking logic here.
    const UInt64 rid = large_run_id & ~(1ULL << 63); // mirrors clearSignBit
    OpenSwathOSWWriter::OSWData rows;
    rows.feature_rows.push_back({
      StringUtils::toStr(rid + 1), StringUtils::toStr(rid), "999", "100.0", "NULL", "100.0",
      "0.0", "90.0", "110.0", "NULL", "NULL"
    });
    writer.writeRows(rows);
  }

  // --- Verify ---
  // Use READWRITE to allow creating temp tables for verification queries
  {
    SqliteConnector conn(temp_file);

    // 1. Check typeof(RUN.ID) is 'integer' (was 'blob' before the fix)
    conn.executeStatement(
      "CREATE TEMP TABLE _typeof_ok AS SELECT 1 FROM RUN WHERE typeof(ID) = 'integer';");
    TEST_EQUAL(conn.countTableRows("_typeof_ok"), 1)

    conn.executeStatement(
      "CREATE TEMP TABLE _feature_typeof_ok AS SELECT 1 FROM FEATURE WHERE typeof(RUN_ID) = 'integer';");
    TEST_EQUAL(conn.countTableRows("_feature_typeof_ok"), 1)

    // 2. Check that JOIN on FEATURE.RUN_ID = RUN.ID returns rows
    conn.executeStatement(
      "CREATE TEMP TABLE _join_result AS "
      "SELECT RUN.ID FROM RUN INNER JOIN FEATURE ON FEATURE.RUN_ID = RUN.ID;");
    TEST_EQUAL(conn.countTableRows("_join_result"), 1)

    // 3. Verify the stored RUN.ID value round-trips correctly
    const UInt64 rid = large_run_id & ~(1ULL << 63);
    conn.executeStatement(
      "CREATE TEMP TABLE _id_match AS SELECT 1 FROM RUN WHERE ID = " + StringUtils::toStr(rid) + ";");
    TEST_EQUAL(conn.countTableRows("_id_match"), 1)
  }

  File::remove(temp_file);
}
END_SECTION

START_SECTION([EXTRA] OSWValue holds numbers in place and text out of line)
{
  typedef OpenSwathOSWWriter::OSWValue OSWValue;
  // a scored feature fills several hundred cells: keep them small
  TEST_EQUAL(sizeof(OSWValue) <= 16, true)

  OSWValue int_value(Int64(42));
  TEST_EQUAL(int_value.type == OSWValue::Type::Int64, true)
  TEST_EQUAL(int_value.asInt(), 42)
  OSWValue double_value(1.5);
  TEST_EQUAL(double_value.type == OSWValue::Type::Double, true)
  TEST_REAL_SIMILAR(double_value.asDouble(), 1.5)
  TEST_EQUAL(OSWValue("NULL").isNull(), true)
  TEST_EQUAL(OSWValue("17").asInt(), 17)

  const std::string long_text = "a text value that does not fit into the small string buffer";
  OSWValue text = OSWValue::text(long_text);
  TEST_EQUAL(text.type == OSWValue::Type::Text, true)
  TEST_EQUAL(text.asText(), long_text)

  // copies own their text
  OSWValue copy(text);
  OSWValue assigned;
  assigned = text;
  text = OSWValue(Int64(1));
  TEST_EQUAL(text.asInt(), 1)
  TEST_EQUAL(copy.asText(), long_text)
  TEST_EQUAL(assigned.asText(), long_text)

  // moving hands the text over and leaves a NULL cell behind
  OSWValue moved(std::move(copy));
  TEST_EQUAL(moved.asText(), long_text)
  TEST_EQUAL(copy.isNull(), true)
  assigned = std::move(moved);
  TEST_EQUAL(assigned.asText(), long_text)
  TEST_EQUAL(moved.isNull(), true)
  assigned = OSWValue::text("short");
  TEST_EQUAL(assigned.asText(), "short")

  // rows of cells can be copied and grown (vector reallocation moves the cells)
  std::vector<OpenSwathOSWWriter::FeaturePrecursorRow> rows;
  for (int i = 0; i < 100; ++i)
  {
    rows.push_back({OSWValue(Int64(i)), OSWValue::text(long_text + StringUtils::toStr(i)), OSWValue(0.5 * i), OSWValue()});
  }
  const std::vector<OpenSwathOSWWriter::FeaturePrecursorRow> rows_copy = rows;
  TEST_EQUAL(rows_copy[99][0].asInt(), 99)
  TEST_EQUAL(rows_copy[99][1].asText(), long_text + "99")
  TEST_REAL_SIMILAR(rows_copy[7][2].asDouble(), 3.5)
  TEST_EQUAL(rows_copy[7][3].isNull(), true)
}
END_SECTION

START_SECTION([EXTRA] prepareLine indexes masserror_ppm by MS2 subordinate position)
{
  std::string temp_file = TempFiles::getTemporaryFile();

  {
    OpenSwathOSWWriter writer(temp_file, false);
    writer.writeHeader();
    writer.addRun(1, "test_input_file.mzML");

    Feature feature;
    feature.ensureUniqueId();
    feature.setRT(10.0);
    feature.setIntensity(100.0);
    feature.setMetaValue("leftWidth", 9.0);
    feature.setMetaValue("rightWidth", 11.0);
    feature.setMetaValue("total_xic", 1000.0);
    feature.setMetaValue("peak_apices_sum", 100.0);
    feature.setMetaValue("masserror_ppm", DataValue(DoubleList{12.5}));

    Feature ms1;
    ms1.setIntensity(10.0);
    ms1.setMetaValue("FeatureLevel", "MS1");
    ms1.setMetaValue("native_id", "Precursor_i0");
    ms1.setMetaValue("peak_apex_int", 5.0);

    Feature ms2;
    ms2.setIntensity(50.0);
    ms2.setMetaValue("FeatureLevel", "MS2");
    ms2.setMetaValue("native_id", 123);
    ms2.setMetaValue("total_xic", 500.0);
    ms2.setMetaValue("peak_apex_int", 25.0);
    ms2.setMetaValue("peak_apex_position", 10.5);
    ms2.setMetaValue("width_at_50", 2.0);

    feature.getSubordinates().push_back(ms1);
    feature.getSubordinates().push_back(ms2);

    FeatureMap output;
    output.push_back(feature);
    writer.writeLines(std::vector<std::string>{writer.prepareLine(OpenSwath::LightCompound(), nullptr, output, "77")});
  }

  {
    SqliteConnector conn(temp_file);
    conn.executeStatement(
      "CREATE TEMP TABLE _masserror_ok AS "
      "SELECT 1 FROM FEATURE_TRANSITION WHERE TRANSITION_ID = 123 AND MASSERROR_PPM = 12.5;");
    TEST_EQUAL(conn.countTableRows("_masserror_ok"), 1)
  }

  File::remove(temp_file);
}
END_SECTION

START_SECTION([EXTRA] prepareLine uses available UIS transition names instead of reported count)
{
  std::string temp_file = TempFiles::getTemporaryFile();

  {
    OpenSwathOSWWriter writer(temp_file, true);
    writer.writeHeader();
    writer.addRun(1, "test_input_file.mzML");

    Feature feature;
    feature.ensureUniqueId();
    feature.setRT(10.0);
    feature.setIntensity(100.0);
    feature.setMetaValue("leftWidth", 9.0);
    feature.setMetaValue("rightWidth", 11.0);
    feature.setMetaValue("total_xic", 1000.0);
    feature.setMetaValue("peak_apices_sum", 100.0);
    feature.setMetaValue("id_target_num_transitions", 2);
    feature.setMetaValue("id_target_transition_names", DataValue(IntList{201}));
    feature.setMetaValue("id_target_area_intensity", DataValue(DoubleList{10.0}));
    feature.setMetaValue("id_target_total_area_intensity", DataValue(DoubleList{20.0}));
    feature.setMetaValue("id_target_apex_intensity", DataValue(DoubleList{30.0}));
    feature.setMetaValue("id_target_peak_apex_position", DataValue(DoubleList{10.5}));
    feature.setMetaValue("id_target_width_at_50", DataValue(DoubleList{2.0}));

    FeatureMap output;
    output.push_back(feature);
    writer.writeLines(std::vector<std::string>{writer.prepareLine(OpenSwath::LightCompound(), nullptr, output, "77")});
  }

  {
    SqliteConnector conn(temp_file);
    conn.executeStatement("CREATE TEMP TABLE _uis_rows AS SELECT 1 FROM FEATURE_TRANSITION;");
    TEST_EQUAL(conn.countTableRows("_uis_rows"), 1)

    conn.executeStatement(
      "CREATE TEMP TABLE _uis_expected AS "
      "SELECT 1 FROM FEATURE_TRANSITION WHERE TRANSITION_ID = 201 AND AREA_INTENSITY = 10.0;");
    TEST_EQUAL(conn.countTableRows("_uis_expected"), 1)
  }

  File::remove(temp_file);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
