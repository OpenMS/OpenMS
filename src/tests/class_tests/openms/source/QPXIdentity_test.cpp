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

#include <OpenMS/FORMAT/QPXIdentity.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

// Every expected value in this file was produced by the reference implementation,
// qpx.core.data.identity (bigbio/qpx, branch dev), NOT by running this code and recording what
// it printed. That is the whole point of the test: OpenMS must derive the identity qpx derives,
// because qpxc re-derives feature_id on conversion while leaving psm.feature_id alone, so ids
// that merely look plausible would turn into dangling cross-references.

START_TEST(QPXIdentity, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((std::string formatFloat(float value)))
{
  // Python's repr() rules, which json.dumps inherits: shortest round-tripping digits, fixed
  // notation unless the decimal point would land at <= -4 or > 16, and ".0" appended so a float
  // never prints as a bare integer.
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(0.0f), "0.0")
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(-0.0f), "-0.0")
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(1.0f), "1.0")
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(-1.0f), "-1.0")
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(123456.0f), "123456.0")

  // A float32 is widened to double before it is rendered, so what gets hashed is the double
  // that Arrow's cast produces -- 0.1f is 0.10000000149011612, not 0.1.
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(0.1f), "0.10000000149011612")
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(549.8573f), "549.8572998046875")
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(1534.9807f), "1534.980712890625")

  // The fixed/exponential switch, probed on both sides of each boundary. Powers of two, so the
  // values are exact in float32 AND double and the test is about the boundary rather than about
  // which decimal a rounding happened to produce.
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(1125899906842624.0f), "1125899906842624.0")   // 2^50, decpt 16: fixed
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(18014398509481984.0f), "1.8014398509481984e+16") // 2^54, decpt 17: exponential
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(0.0001220703125f), "0.0001220703125")         // 2^-13, decpt -3: fixed
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(0.00006103515625f), "6.103515625e-05")        // 2^-14, decpt -4: exponential
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(1e-4f), "9.999999747378752e-05")

  // Subnormal and maximum float32; the exponent is always signed and at least two digits.
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(1e-45f), "1.401298464324817e-45")
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(3.4028235e38f), "3.4028234663852886e+38")

  // json.dumps writes these bare, not quoted. Unreachable from a validated QPX table, but the
  // encoder must not produce something unparseable if it ever sees one.
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(std::numeric_limits<float>::infinity()), "Infinity")
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(-std::numeric_limits<float>::infinity()), "-Infinity")
  TEST_STRING_EQUAL(QPXIdentity::formatFloat(std::numeric_limits<float>::quiet_NaN()), "NaN")
}
END_SECTION

START_SECTION((std::string canonical(const std::vector<Value>& values, const std::vector<size_t>& unordered_list_indices)))
{
  using V = QPXIdentity::Value;

  TEST_STRING_EQUAL(QPXIdentity::canonical({V{std::string("BSA1_F1")}, V{std::vector<Int32>{2311}},
                                            V{std::string("PEPTIDER")}, V{Int64(2)}}),
                    "[\"BSA1_F1\",[2311],\"PEPTIDER\",2]")

  // Compact separators: no space after ',' or ':'.
  TEST_STRING_EQUAL(QPXIdentity::canonical({V{Int64(1)}, V{Int64(2)}}), "[1,2]")

  // A null column value is JSON null, not an empty string -- the two must not collide.
  TEST_STRING_EQUAL(QPXIdentity::canonical({V{QPXIdentity::Null{}}}), "[null]")
  TEST_NOT_EQUAL(QPXIdentity::canonical({V{QPXIdentity::Null{}}}),
                 QPXIdentity::canonical({V{std::string("")}}))

  // ensure_ascii: non-ASCII is escaped, in lower-case hex.
  // NB the "" after \xa4: a C++ hex escape is greedy, so \xa4cker would be read as one number.
  TEST_STRING_EQUAL(QPXIdentity::canonical({V{std::string("B\xc3\xa4""cker")}}), "[\"B\\u00e4cker\"]")
  // Outside the BMP, as a surrogate pair (U+1F600).
  TEST_STRING_EQUAL(QPXIdentity::canonical({V{std::string("\xf0\x9f\x98\x80")}}), "[\"\\ud83d\\ude00\"]")
  // Escapes required by JSON itself.
  TEST_STRING_EQUAL(QPXIdentity::canonical({V{std::string("a\"b\\c\nd")}}), "[\"a\\\"b\\\\c\\nd\"]")
  // DEL is not printable ASCII and is escaped; '/' is not escaped.
  TEST_STRING_EQUAL(QPXIdentity::canonical({V{std::string("a/b\x7f")}}), "[\"a/b\\u007f\"]")

  // JSON quoting is what stops distinct composites from aliasing onto the same bytes: a plain
  // delimiter join would make these two indistinguishable.
  TEST_NOT_EQUAL(QPXIdentity::canonical({V{std::string("a,b")}, V{std::string("c")}}),
                 QPXIdentity::canonical({V{std::string("a")}, V{std::string("b,c")}}))

  // Set-valued fields are sorted; ordered ones are not.
  const std::vector<std::string> runs_forward{"BSA1_F1", "BSA1_F2"};
  const std::vector<std::string> runs_reversed{"BSA1_F2", "BSA1_F1"};
  TEST_STRING_EQUAL(QPXIdentity::canonical({V{runs_reversed}}, {0}), "[[\"BSA1_F1\",\"BSA1_F2\"]]")
  TEST_STRING_EQUAL(QPXIdentity::canonical({V{runs_forward}}, {0}),
                    QPXIdentity::canonical({V{runs_reversed}}, {0}))
  TEST_NOT_EQUAL(QPXIdentity::canonical({V{runs_forward}}), QPXIdentity::canonical({V{runs_reversed}}))
  TEST_STRING_EQUAL(QPXIdentity::canonical({V{std::vector<Int32>{3, 1, 2}}}), "[[3,1,2]]")
}
END_SECTION

START_SECTION((Int64 deriveId(const std::vector<Value>& values, const std::vector<size_t>& unordered_list_indices)))
{
  using V = QPXIdentity::Value;

  // qpx: derive_id(["BSA1_F1", [2311], "PEPTIDER", 2])
  TEST_EQUAL(QPXIdentity::deriveId({V{std::string("BSA1_F1")}, V{std::vector<Int32>{2311}},
                                    V{std::string("PEPTIDER")}, V{Int64(2)}}),
             Int64(249030645178071364))

  // Negative values are normal: the digest is read as a signed big-endian int64.
  TEST_EQUAL(QPXIdentity::deriveId({V{std::string("B\xc3\xa4""cker/prote\xc3\xafn")},
                                    V{std::vector<Int32>{1, 2, 3}},
                                    V{std::string("PEP[Phospho]TIDE")}, V{Int64(3)}}),
             Int64(-3403684637348881565))

  // Composites whose text would collide under a delimiter join stay distinct here.
  TEST_EQUAL(QPXIdentity::deriveId({V{std::string("run,with]brackets")}, V{std::vector<Int32>{}},
                                    V{std::string("")}, V{Int64(0)}}),
             Int64(260970448792327776))
}
END_SECTION

START_SECTION((Int64 featureId(const std::string& run_file_name, const std::string& peptidoform, Int64 charge, std::optional<float> rt, const std::vector<Int32>& scan, float observed_mz)))
{
  // qpx: derive_id(["BSA1_F1", "PEPTIDER", 2, 1534.980712890625, [2311], 549.8572998046875])
  TEST_EQUAL(QPXIdentity::featureId("BSA1_F1", "PEPTIDER", 2, 1534.9807f, {2311}, 549.8573f),
             Int64(-9185367007480417157))

  // An unidentified feature row: empty peptidoform, empty scan. Only rt and observed_mz keep it
  // apart from its neighbours, which is why the composite carries floats at all.
  TEST_EQUAL(QPXIdentity::featureId("BSA1_F1", "", 2, 1534.9807f, {}, 549.8573f),
             Int64(6679072350175538872))

  // Each component actually participates.
  const Int64 base = QPXIdentity::featureId("BSA1_F1", "PEPTIDER", 2, 1534.9807f, {2311}, 549.8573f);
  TEST_NOT_EQUAL(QPXIdentity::featureId("BSA1_F2", "PEPTIDER", 2, 1534.9807f, {2311}, 549.8573f), base)
  TEST_NOT_EQUAL(QPXIdentity::featureId("BSA1_F1", "PEPTIDEK", 2, 1534.9807f, {2311}, 549.8573f), base)
  TEST_NOT_EQUAL(QPXIdentity::featureId("BSA1_F1", "PEPTIDER", 3, 1534.9807f, {2311}, 549.8573f), base)
  TEST_NOT_EQUAL(QPXIdentity::featureId("BSA1_F1", "PEPTIDER", 2, 1534.9808f, {2311}, 549.8573f), base)
  TEST_NOT_EQUAL(QPXIdentity::featureId("BSA1_F1", "PEPTIDER", 2, 1534.9807f, {2312}, 549.8573f), base)
  TEST_NOT_EQUAL(QPXIdentity::featureId("BSA1_F1", "PEPTIDER", 2, 1534.9807f, {2311}, 549.8574f), base)

  // A null rt is a key value of its own, distinct from any number.
  TEST_NOT_EQUAL(QPXIdentity::featureId("BSA1_F1", "PEPTIDER", 2, std::nullopt, {2311}, 549.8573f), base)

  // The identity is a function of the PERSISTED float32, so a double that narrows to the same
  // float32 must give the same id -- otherwise ids would depend on the caller's arithmetic.
  TEST_EQUAL(QPXIdentity::featureId("BSA1_F1", "PEPTIDER", 2,
                                    static_cast<float>(1534.98071289062500001), {2311}, 549.8573f),
             base)
}
END_SECTION

START_SECTION((Int64 psmId(const std::string& run_file_name, const std::vector<Int32>& scan, const std::string& peptidoform, Int64 charge)))
{
  // qpx: derive_id(["BSA1_F1", [2311], "PEPTIDER", 2]) -- the psm composite is the same four
  // values in the same order, so it must agree with the deriveId() vector above.
  TEST_EQUAL(QPXIdentity::psmId("BSA1_F1", {2311}, "PEPTIDER", 2), Int64(249030645178071364))

  // Multi-component scans keep their order: a psm identity is not a set.
  TEST_NOT_EQUAL(QPXIdentity::psmId("BSA1_F1", {1, 2}, "PEPTIDER", 2),
                 QPXIdentity::psmId("BSA1_F1", {2, 1}, "PEPTIDER", 2))

  // No floats in this composite, which is what makes psm ids stable across platforms.
  TEST_NOT_EQUAL(QPXIdentity::psmId("BSA1_F1", {2311}, "PEPTIDER", 3),
                 QPXIdentity::psmId("BSA1_F1", {2311}, "PEPTIDER", 2))
}
END_SECTION

START_SECTION((Int64 pgId(const std::string& anchor_protein, const std::vector<std::string>& grouped_runs, const std::optional<std::string>& label)))
{
  // qpx: derive_id(["P02769", ["BSA1_F1", "BSA1_F2"], "LFQ"], unordered_list_indices=(1,))
  TEST_EQUAL(QPXIdentity::pgId("P02769", {"BSA1_F1", "BSA1_F2"}, std::string("LFQ")),
             Int64(-3984318991531137877))

  // grouped_runs is a set: the order the design happens to list the fractions in must not
  // change the group's identity.
  TEST_EQUAL(QPXIdentity::pgId("P02769", {"BSA1_F2", "BSA1_F1"}, std::string("LFQ")),
             Int64(-3984318991531137877))

  // An identification-only group carries no quantity; its null label is part of the key.
  TEST_EQUAL(QPXIdentity::pgId("P02769", {"BSA1_F1"}, std::nullopt), Int64(6696460256469222492))
  TEST_NOT_EQUAL(QPXIdentity::pgId("P02769", {"BSA1_F1"}, std::string("")),
                 QPXIdentity::pgId("P02769", {"BSA1_F1"}, std::nullopt))

  // Different membership, different group -- even when the runs are a subset.
  TEST_NOT_EQUAL(QPXIdentity::pgId("P02769", {"BSA1_F1"}, std::string("LFQ")),
                 QPXIdentity::pgId("P02769", {"BSA1_F1", "BSA1_F2"}, std::string("LFQ")))
}
END_SECTION

START_SECTION([EXTRA] identity_composite footer declarations)
{
  // These strings go into the Parquet footer, and qpxc reads them back to re-derive the ids.
  // They must name the columns in the order the values are hashed in.
  TEST_STRING_EQUAL(QPXIdentity::FEATURE_COMPOSITE,
                    "run_file_name,peptidoform,charge,rt,scan,observed_mz")
  TEST_STRING_EQUAL(QPXIdentity::PSM_COMPOSITE, "run_file_name,scan,peptidoform,charge")
  TEST_STRING_EQUAL(QPXIdentity::PG_COMPOSITE, "anchor_protein,grouped_runs,label")
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
