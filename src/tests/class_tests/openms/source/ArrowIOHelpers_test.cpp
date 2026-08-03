// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
///////////////////////////

#include <OpenMS/KERNEL/ConsensusMap.h>

#include <arrow/api.h>

#include <set>

using namespace OpenMS;
using namespace std;

START_TEST(ArrowIOHelpers, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((std::string generateUuidV4()))
{
  std::string u = ArrowIOHelpers::generateUuidV4();

  // canonical 8-4-4-4-12 form with hyphens
  TEST_EQUAL(u.size(), 36)
  TEST_EQUAL(u[8], '-')
  TEST_EQUAL(u[13], '-')
  TEST_EQUAL(u[18], '-')
  TEST_EQUAL(u[23], '-')

  // version nibble must be '4' (UUID v4)
  TEST_EQUAL(u[14], '4')

  // variant nibble encodes the 10xx bits -> one of 8, 9, a, b
  const char var = u[19];
  TEST_EQUAL(var == '8' || var == '9' || var == 'a' || var == 'b' ||
             var == 'A' || var == 'B', true)

  // every non-hyphen character is a hex digit
  bool all_hex = true;
  for (size_t i = 0; i < u.size(); ++i)
  {
    if (i == 8 || i == 13 || i == 18 || i == 23) continue;
    const char c = u[i];
    const bool hex = (c >= '0' && c <= '9') || (c >= 'a' && c <= 'f') || (c >= 'A' && c <= 'F');
    if (!hex) { all_hex = false; break; }
  }
  TEST_EQUAL(all_hex, true)

  // uniqueness across many draws (122 random bits -> collisions negligible)
  std::set<std::string> seen;
  for (int i = 0; i < 1000; ++i) seen.insert(ArrowIOHelpers::generateUuidV4());
  TEST_EQUAL(seen.size(), 1000)
}
END_SECTION

START_SECTION([EXTRA] typed value accessors getColumn and isNull)
{
  // Build a 2-row table: row 0 carries a value, row 1 is null in every column.
  std::shared_ptr<arrow::Array> s_arr, d_arr, f_arr, i32_arr, i64_arr, b_arr;
  {
    arrow::StringBuilder b;
    TEST_EQUAL(b.Append("hello").ok(), true)
    TEST_EQUAL(b.AppendNull().ok(), true)
    TEST_EQUAL(b.Finish(&s_arr).ok(), true)
  }
  {
    arrow::DoubleBuilder b;
    TEST_EQUAL(b.Append(2.5).ok(), true)
    TEST_EQUAL(b.AppendNull().ok(), true)
    TEST_EQUAL(b.Finish(&d_arr).ok(), true)
  }
  {
    arrow::FloatBuilder b;
    TEST_EQUAL(b.Append(1.25f).ok(), true)
    TEST_EQUAL(b.AppendNull().ok(), true)
    TEST_EQUAL(b.Finish(&f_arr).ok(), true)
  }
  {
    arrow::Int32Builder b;
    TEST_EQUAL(b.Append(42).ok(), true)
    TEST_EQUAL(b.AppendNull().ok(), true)
    TEST_EQUAL(b.Finish(&i32_arr).ok(), true)
  }
  {
    arrow::Int64Builder b;
    TEST_EQUAL(b.Append(static_cast<int64_t>(9000000000LL)).ok(), true)
    TEST_EQUAL(b.AppendNull().ok(), true)
    TEST_EQUAL(b.Finish(&i64_arr).ok(), true)
  }
  {
    arrow::BooleanBuilder b;
    TEST_EQUAL(b.Append(true).ok(), true)
    TEST_EQUAL(b.AppendNull().ok(), true)
    TEST_EQUAL(b.Finish(&b_arr).ok(), true)
  }

  // row 0: values are read back
  TEST_EQUAL(ArrowIOHelpers::getStringValue(s_arr, 0), "hello")
  TEST_REAL_SIMILAR(ArrowIOHelpers::getDoubleValue(d_arr, 0), 2.5)
  TEST_REAL_SIMILAR(ArrowIOHelpers::getFloatValue(f_arr, 0), 1.25f)
  TEST_EQUAL(ArrowIOHelpers::getInt32Value(i32_arr, 0), 42)
  TEST_EQUAL(ArrowIOHelpers::getInt64Value(i64_arr, 0), static_cast<int64_t>(9000000000LL))
  TEST_EQUAL(ArrowIOHelpers::getBoolValue(b_arr, 0), true)

  // row 1: null -> default value / empty string; isNull reports true
  TEST_EQUAL(ArrowIOHelpers::getStringValue(s_arr, 1), "")
  TEST_REAL_SIMILAR(ArrowIOHelpers::getDoubleValue(d_arr, 1, -1.0), -1.0)
  TEST_EQUAL(ArrowIOHelpers::getInt64Value(i64_arr, 1, -7), static_cast<int64_t>(-7))
  TEST_EQUAL(ArrowIOHelpers::getBoolValue(b_arr, 1, true), true)  // default returned
  TEST_EQUAL(ArrowIOHelpers::isNull(s_arr, 1), true)
  TEST_EQUAL(ArrowIOHelpers::isNull(s_arr, 0), false)

  // a null array pointer is treated as null everywhere
  std::shared_ptr<arrow::Array> none;
  TEST_EQUAL(ArrowIOHelpers::isNull(none, 0), true)
  TEST_EQUAL(ArrowIOHelpers::getStringValue(none, 0), "")
  TEST_EQUAL(ArrowIOHelpers::getInt32Value(none, 0, 99), 99)

  // getColumn retrieves a column by name; a missing column (required=false) -> null
  auto schema = arrow::schema({arrow::field("str", arrow::utf8()),
                               arrow::field("i64", arrow::int64())});
  auto table = arrow::Table::Make(schema, {s_arr, i64_arr});

  auto col_str = ArrowIOHelpers::getColumn(table, "str");
  TEST_EQUAL(col_str != nullptr, true)
  TEST_EQUAL(ArrowIOHelpers::getStringValue(col_str, 0), "hello")

  auto col_i64 = ArrowIOHelpers::getColumn(table, "i64");
  TEST_EQUAL(col_i64 != nullptr, true)
  TEST_EQUAL(ArrowIOHelpers::getInt64Value(col_i64, 0), static_cast<int64_t>(9000000000LL))

  auto missing = ArrowIOHelpers::getColumn(table, "does_not_exist", false);
  TEST_EQUAL(missing == nullptr, true)
}
END_SECTION

START_SECTION((std::shared_ptr<const arrow::KeyValueMetadata> qpxFileMetadata(const std::string&, const ParquetWriteConfig&, const std::map<std::string, std::string>&)))
{
  // Default config is ZSTD -> the full documented key set is present.
  auto md = ArrowIOHelpers::qpxFileMetadata("pg_file");
  TEST_NOT_EQUAL(md, nullptr)

  auto value_of = [&md](const std::string& key) -> std::string
  {
    auto r = md->Get(key);
    return r.ok() ? r.ValueOrDie() : std::string();
  };

  TEST_STRING_EQUAL(value_of("qpx_version"), "1.1")
  TEST_STRING_EQUAL(value_of("file_type"), "pg_file")
  TEST_STRING_EQUAL(value_of("creator"), "OpenMS")
  TEST_STRING_EQUAL(value_of("compression_format"), "zstd")
  // software_provider carries the version; creator stays the bare org/tool name
  TEST_TRUE(StringUtils::hasPrefix(value_of("software_provider"), "OpenMS "))
  TEST_EQUAL(value_of("uuid").size(), 36)
  TEST_TRUE(!value_of("creation_date").empty())

  // Each call is a distinct file identity — reusing one call per file is the caller's job.
  auto md2 = ArrowIOHelpers::qpxFileMetadata("pg_file");
  TEST_NOT_EQUAL(md2, nullptr)
  TEST_NOT_EQUAL(value_of("uuid"), md2->Get("uuid").ValueOrDie())

  // Extra keys are appended verbatim.
  auto md3 = ArrowIOHelpers::qpxFileMetadata("psm_file", ParquetWriteConfig{}, {{"scan_format", "index"}});
  TEST_NOT_EQUAL(md3, nullptr)
  TEST_STRING_EQUAL(md3->Get("scan_format").ValueOrDie(), "index")

  // compression_format vocabulary is zstd|snappy|gzip|lzo|none.
  ParquetWriteConfig cfg;
  cfg.compression = ParquetWriteConfig::Compression::SNAPPY;
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxFileMetadata("psm_file", cfg)->Get("compression_format").ValueOrDie(), "snappy")
  cfg.compression = ParquetWriteConfig::Compression::GZIP;
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxFileMetadata("psm_file", cfg)->Get("compression_format").ValueOrDie(), "gzip")
  cfg.compression = ParquetWriteConfig::Compression::NONE;
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxFileMetadata("psm_file", cfg)->Get("compression_format").ValueOrDie(), "none")

  // LZ4 has no QPX token: refuse rather than emit an out-of-vocabulary value.
  cfg.compression = ParquetWriteConfig::Compression::LZ4;
  TEST_EQUAL(ArrowIOHelpers::qpxFileMetadata("psm_file", cfg), nullptr)
}
END_SECTION

START_SECTION((bool qpxWarnOnRunNameCollisions(const std::string&, const std::vector<std::string>&)))
{
  // Warn-only, matching the PSM exporter: same-named files in different directories are a
  // legitimate layout, so an ambiguous key is reported rather than refused.
  TEST_TRUE(ArrowIOHelpers::qpxWarnOnRunNameCollisions("t", {"/a/run1.mzML", "/a/run2.mzML"}))
  TEST_TRUE(ArrowIOHelpers::qpxWarnOnRunNameCollisions("t", {"/a/run1.mzML", "/a/run1.mzML"}))  // one run, twice
  TEST_TRUE(ArrowIOHelpers::qpxWarnOnRunNameCollisions("t", {}))
  TEST_TRUE(ArrowIOHelpers::qpxWarnOnRunNameCollisions("t", {"", "/a/run1.mzML"}))

  // Distinct paths collapsing onto one run_file_name: reported (false), not fatal.
  TEST_FALSE(ArrowIOHelpers::qpxWarnOnRunNameCollisions("t", {"/frac_a/run1.mzML", "/frac_b/run1.mzML"}))
  TEST_FALSE(ArrowIOHelpers::qpxWarnOnRunNameCollisions("t", {"/a/run1.mzML", "/a/run1.raw"}))
}
END_SECTION

START_SECTION((std::string qpxIntensityLabel(const std::string&, const std::string&)))
{
  // Label-free: ProteomicsLFQ stamps "label-free" on every header; an unset label is the same.
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("label-free", ""), "LFQ")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("", ""), "LFQ")

  // Isobaric: IsobaricChannelExtractor builds the label as "<methodname>_<channelname>".
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("tmt6plex_126", "126"), "TMT126")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("tmt10plex_126", "126"), "TMT126")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("tmt10plex_127N", "127N"), "TMT127N")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("tmt11plex_131C", "131C"), "TMT131C")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("tmt16plex_134N", "134N"), "TMT134N")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("tmt18plex_135N", "135N"), "TMT135N")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("itraq4plex_114", "114"), "ITRAQ114")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("itraq8plex_121", "121"), "ITRAQ121")

  // Older and synthetic maps may encode the complete identity only in the label.
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("tmt6plex_126", ""), "TMT126")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("itraq4plex_114", ""), "ITRAQ114")

  // TMT10-plex channel 10 is "131" in OpenMS' naming. qpx's own converter map is
  // 11-plex-indexed and calls the 10th channel "TMT131N"; OpenMS' name is authoritative.
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("tmt10plex_131", "131"), "TMT131")
  // ... and 11-plex really does have both 131N and 131C.
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("tmt11plex_131N", "131N"), "TMT131N")

  // Scalar SILAC normalization uses the modification's mass class. Whole-map normalization below
  // additionally sees the plex role, which matters for an Arg6 two-plex.
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("no_label", ""), "SILAC light")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("Arg6", ""), "SILAC medium")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("Lys8Arg10", ""), "SILAC heavy")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("Dimethyl4", ""), "DIMETHYL4")

  // An unstandardized multiplex token cannot be written into a QPX join key.
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("label 0", ""), "")

  // A channel whose method cannot be identified must NOT be guessed — it is a join key.
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("mystery_126", "126"), "")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxIntensityLabel("no_separator", "126"), "")
}
END_SECTION

START_SECTION((std::map<std::uint64_t, std::string> qpxIntensityLabels(const ConsensusMap&)))
{
  ConsensusMap two_plex;
  two_plex.setExperimentType("labeled_MS1");
  ConsensusMap::ColumnHeader light;
  light.filename = "silac.mzML";
  light.label = "no_label";
  light.setMetaValue("channel_id", 0);
  ConsensusMap::ColumnHeader arg6;
  arg6.filename = "silac.mzML";
  arg6.label = "Arg6";
  arg6.setMetaValue("channel_id", 1);
  two_plex.setColumnHeaders({{0, light}, {1, arg6}});

  auto labels = ArrowIOHelpers::qpxIntensityLabels(two_plex);
  TEST_STRING_EQUAL(labels.at(0), "SILAC light")
  // Arg6 is medium by mass class, but it is the non-light (heavy-role) channel of a two-plex.
  TEST_STRING_EQUAL(labels.at(1), "SILAC heavy")

  ConsensusMap three_plex;
  three_plex.setExperimentType("labeled_MS1");
  ConsensusMap::ColumnHeader heavy;
  heavy.filename = "silac3.mzML";
  heavy.label = "Lys8Arg10";
  heavy.setMetaValue("channel_id", 2);
  ConsensusMap::ColumnHeader medium;
  medium.filename = "silac3.mzML";
  medium.label = "Lys4Arg6";
  medium.setMetaValue("channel_id", 1);
  light.filename = "silac3.mzML";
  three_plex.setColumnHeaders({{3, heavy}, {7, light}, {9, medium}});

  labels = ArrowIOHelpers::qpxIntensityLabels(three_plex);
  TEST_STRING_EQUAL(labels.at(3), "SILAC heavy")
  TEST_STRING_EQUAL(labels.at(7), "SILAC light")
  TEST_STRING_EQUAL(labels.at(9), "SILAC medium")

  ConsensusMap unsupported;
  unsupported.setExperimentType("labeled_MS1");
  heavy.filename = "unsupported_silac.mzML";
  unsupported.setColumnHeaders({{0, heavy}});
  labels = ArrowIOHelpers::qpxIntensityLabels(unsupported);
  TEST_STRING_EQUAL(labels.at(0), "")

  ConsensusMap unknown_role;
  unknown_role.setExperimentType("labeled_MS1");
  light.filename = "unknown_role.mzML";
  ConsensusMap::ColumnHeader unknown = heavy;
  unknown.filename = "unknown_role.mzML";
  unknown.label = "SILAC mystery";
  unknown_role.setColumnHeaders({{0, light}, {1, unknown}});
  labels = ArrowIOHelpers::qpxIntensityLabels(unknown_role);
  TEST_STRING_EQUAL(labels.at(0), "")
  TEST_STRING_EQUAL(labels.at(1), "")

  // A two-plex needs exactly one light channel.
  ConsensusMap no_light;
  no_light.setExperimentType("labeled_MS1");
  medium.filename = "no_light.mzML";
  heavy.filename = "no_light.mzML";
  no_light.setColumnHeaders({{0, medium}, {1, heavy}});
  labels = ArrowIOHelpers::qpxIntensityLabels(no_light);
  TEST_STRING_EQUAL(labels.at(0), "")
  TEST_STRING_EQUAL(labels.at(1), "")

  // A three-plex needs one channel per role.
  ConsensusMap duplicate_role;
  duplicate_role.setExperimentType("labeled_MS1");
  light.filename = "duplicate_role.mzML";
  medium.filename = "duplicate_role.mzML";
  ConsensusMap::ColumnHeader medium_again = medium;
  duplicate_role.setColumnHeaders({{0, light}, {1, medium}, {2, medium_again}});
  labels = ArrowIOHelpers::qpxIntensityLabels(duplicate_role);
  TEST_STRING_EQUAL(labels.at(0), "")
  TEST_STRING_EQUAL(labels.at(1), "")
  TEST_STRING_EQUAL(labels.at(2), "")

  ConsensusMap label_free;
  label_free.setExperimentType("label-free");
  ConsensusMap::ColumnHeader lfq;
  lfq.filename = "lfq.mzML";
  lfq.label = "label-free";
  label_free.setColumnHeaders({{0, lfq}});
  labels = ArrowIOHelpers::qpxIntensityLabels(label_free);
  TEST_STRING_EQUAL(labels.at(0), "LFQ")

  // Ordinary words containing "lys" must not route the whole source through SILAC validation.
  ConsensusMap ordinary_labels;
  ordinary_labels.setExperimentType("label-free");
  lfq.filename = "ordinary.mzML";
  ConsensusMap::ColumnHeader lysate = lfq;
  lysate.label = "Lysate_A";
  ConsensusMap::ColumnHeader catalyst = lfq;
  catalyst.label = "catalyst";
  ordinary_labels.setColumnHeaders({{0, lfq}, {1, lysate}, {2, catalyst}});
  labels = ArrowIOHelpers::qpxIntensityLabels(ordinary_labels);
  TEST_STRING_EQUAL(labels.at(0), "LFQ")
  TEST_STRING_EQUAL(labels.at(1), "")
  TEST_STRING_EQUAL(labels.at(2), "")
}
END_SECTION

START_SECTION((bool qpxIsCanonicalIntensityLabel(const std::string&)))
{
  TEST_TRUE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("LFQ"))
  TEST_TRUE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("SILAC light"))
  TEST_TRUE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("SILAC medium"))
  TEST_TRUE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("SILAC heavy"))
  TEST_TRUE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("TMT131"))
  TEST_TRUE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("TMT131N"))
  TEST_TRUE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("ITRAQ121"))
  TEST_TRUE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("DIMETHYL8"))

  TEST_FALSE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("Arg10"))
  TEST_FALSE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("silac heavy"))
  TEST_FALSE(ArrowIOHelpers::qpxIsCanonicalIntensityLabel("TMT999"))
}
END_SECTION

START_SECTION((std::string qpxRunFileName(const std::string&)))
{
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxRunFileName("/data/proj/S1_Frontal_1.mzML"), "S1_Frontal_1")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxRunFileName("S1_Frontal_1.mzML"), "S1_Frontal_1")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxRunFileName("/data/run.raw"), "run")
  // Bruker .d directories reduce like any other path
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxRunFileName("/data/run.d"), "run")
  // Already in QPX form -> unchanged
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxRunFileName("BSA1_F1"), "BSA1_F1")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxRunFileName(""), "")
}
END_SECTION

START_SECTION((std::string qpxScanFormat(const std::string&)))
{
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat("controllerType=0 controllerNumber=1 scan=1234"), "scan")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat("scan=42"), "scan")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat("scanId=42"), "scan")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat("spectrum=7"), "scan")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat("index=7"), "index")
  // Not a recognized native ID -> no evidence either way
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat("some_opaque_identifier"), "")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat(""), "")
}
END_SECTION

START_SECTION((std::string qpxScanFormat(const std::vector<std::string>&)))
{
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat(std::vector<std::string>{"scan=1", "scan=2"}), "scan")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat(std::vector<std::string>{"index=1", "index=2"}), "index")

  // Unrecognized entries carry no evidence and must not veto a consistent set.
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat(std::vector<std::string>{"opaque", "scan=2", ""}), "scan")

  // Genuinely mixed conventions: neither token describes the file, so omit rather than guess.
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat(std::vector<std::string>{"scan=1", "index=2"}), "")

  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat(std::vector<std::string>{}), "")
  TEST_STRING_EQUAL(ArrowIOHelpers::qpxScanFormat(std::vector<std::string>{"opaque", "also_opaque"}), "")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
