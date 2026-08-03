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

#include <OpenMS/METADATA/IdentifierMSRunMapper.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CONCEPT/Constants.h>

///////////////////////////

START_TEST(IdentifierMSRunMapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// Helper: build a small set of ProteinIdentifications with known
// identifiers and primary MS-run-paths.
auto makeProtIds = []() -> std::vector<ProteinIdentification>
{
  ProteinIdentification p1, p2;
  p1.setIdentifier("ID_A");
  p1.setPrimaryMSRunPath(StringList{"/data/runA.mzML"});
  p2.setIdentifier("ID_B");
  p2.setPrimaryMSRunPath(StringList{"/data/runB1.mzML", "/data/runB2.mzML"});
  return std::vector<ProteinIdentification>{p1, p2};
};

IdentifierMSRunMapper* ptr = nullptr;
IdentifierMSRunMapper* null_ptr = nullptr;

START_SECTION((IdentifierMSRunMapper()))
{
  ptr = new IdentifierMSRunMapper();
  TEST_NOT_EQUAL(ptr, null_ptr)
  TEST_EQUAL(ptr->empty(), true)
  TEST_EQUAL(ptr->size(), 0)
}
END_SECTION

START_SECTION((~IdentifierMSRunMapper()))
{
  delete ptr;
}
END_SECTION

START_SECTION((explicit IdentifierMSRunMapper(const std::vector<ProteinIdentification>& prot_ids)))
{
  IdentifierMSRunMapper m(makeProtIds());
  TEST_EQUAL(m.empty(), false)
  TEST_EQUAL(m.size(), 2)
}
END_SECTION

START_SECTION((void create(const std::vector<ProteinIdentification>& prot_ids)))
{
  IdentifierMSRunMapper m;
  m.create(makeProtIds());
  TEST_EQUAL(m.size(), 2)
  TEST_EQUAL(m.hasIdentifier("ID_A"), true)

  // create() clears previous content
  ProteinIdentification p;
  p.setIdentifier("ONLY");
  p.setPrimaryMSRunPath(StringList{"/data/only.mzML"});
  m.create(std::vector<ProteinIdentification>{p});
  TEST_EQUAL(m.size(), 1)
  TEST_EQUAL(m.hasIdentifier("ID_A"), false)
  TEST_EQUAL(m.hasIdentifier("ONLY"), true)

  // Two protein IDs with the same ms-run-path must be rejected
  ProteinIdentification d1, d2;
  d1.setIdentifier("D1");
  d1.setPrimaryMSRunPath(StringList{"/dup.mzML"});
  d2.setIdentifier("D2");
  d2.setPrimaryMSRunPath(StringList{"/dup.mzML"});
  TEST_EXCEPTION(Exception::InvalidValue, m.create(std::vector<ProteinIdentification>{d1, d2}))
}
END_SECTION

START_SECTION((bool empty() const))
{
  IdentifierMSRunMapper m;
  TEST_EQUAL(m.empty(), true)
  m.create(makeProtIds());
  TEST_EQUAL(m.empty(), false)
}
END_SECTION

START_SECTION((Size size() const))
{
  IdentifierMSRunMapper m(makeProtIds());
  TEST_EQUAL(m.size(), 2)
}
END_SECTION

START_SECTION((bool hasIdentifier(const std::string& identifier) const))
{
  IdentifierMSRunMapper m(makeProtIds());
  TEST_EQUAL(m.hasIdentifier("ID_A"), true)
  TEST_EQUAL(m.hasIdentifier("ID_B"), true)
  TEST_EQUAL(m.hasIdentifier("ID_X"), false)
}
END_SECTION

START_SECTION((const StringList& getMSRunPaths(const std::string& identifier) const))
{
  IdentifierMSRunMapper m(makeProtIds());
  const StringList& a = m.getMSRunPaths("ID_A");
  TEST_EQUAL(a.size(), 1)
  TEST_STRING_EQUAL(a[0], "/data/runA.mzML")

  const StringList& b = m.getMSRunPaths("ID_B");
  TEST_EQUAL(b.size(), 2)
  TEST_STRING_EQUAL(b[0], "/data/runB1.mzML")
  TEST_STRING_EQUAL(b[1], "/data/runB2.mzML")

  // Unknown identifier yields an empty list (no exception)
  const StringList& none = m.getMSRunPaths("nope");
  TEST_EQUAL(none.empty(), true)
}
END_SECTION

START_SECTION((std::vector<std::string> getIdentifiers() const))
{
  IdentifierMSRunMapper m(makeProtIds());
  std::vector<std::string> ids = m.getIdentifiers();
  TEST_EQUAL(ids.size(), 2)
  // std::map ordering -> sorted by key
  TEST_STRING_EQUAL(ids[0], "ID_A")
  TEST_STRING_EQUAL(ids[1], "ID_B")
}
END_SECTION

START_SECTION((const std::string& getIdentifier(const StringList& ms_run_paths) const))
{
  IdentifierMSRunMapper m(makeProtIds());
  TEST_STRING_EQUAL(m.getIdentifier(StringList{"/data/runA.mzML"}), "ID_A")
  TEST_STRING_EQUAL(m.getIdentifier(StringList{"/data/runB1.mzML", "/data/runB2.mzML"}), "ID_B")

  // Lookup is on the full path list: a partial match does not resolve
  TEST_EXCEPTION(Exception::ElementNotFound, m.getIdentifier(StringList{"/data/runB1.mzML"}))
  TEST_EXCEPTION(Exception::ElementNotFound, m.getIdentifier(StringList{"/missing.mzML"}))
}
END_SECTION

START_SECTION((bool hasRunPath(const StringList& ms_run_paths) const))
{
  IdentifierMSRunMapper m(makeProtIds());
  TEST_EQUAL(m.hasRunPath(StringList{"/data/runA.mzML"}), true)
  TEST_EQUAL(m.hasRunPath(StringList{"/data/runB1.mzML", "/data/runB2.mzML"}), true)
  TEST_EQUAL(m.hasRunPath(StringList{"/data/runB1.mzML"}), false)
  TEST_EQUAL(m.hasRunPath(StringList{"/missing.mzML"}), false)
}
END_SECTION

START_SECTION((bool tryGetIdentifier(const StringList& ms_run_paths, std::string& identifier) const))
{
  IdentifierMSRunMapper m(makeProtIds());
  std::string id;
  TEST_EQUAL(m.tryGetIdentifier(StringList{"/data/runA.mzML"}, id), true)
  TEST_STRING_EQUAL(id, "ID_A")

  std::string id2 = "unchanged";
  TEST_EQUAL(m.tryGetIdentifier(StringList{"/missing.mzML"}, id2), false)
  TEST_STRING_EQUAL(id2, "unchanged")
}
END_SECTION

START_SECTION((std::string getPrimaryMSRunPath(const PeptideIdentification& pepid) const))
{
  IdentifierMSRunMapper m(makeProtIds());

  // No merge index -> default to first path
  PeptideIdentification pep;
  pep.setIdentifier("ID_B");
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(pep), "/data/runB1.mzML")

  // Valid merge index -> corresponding path
  pep.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, 1);
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(pep), "/data/runB2.mzML")

  // Out-of-range merge index -> empty string
  pep.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, 5);
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(pep), "")

  // Unknown identifier -> empty string
  PeptideIdentification pep_unknown;
  pep_unknown.setIdentifier("DOES_NOT_EXIST");
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(pep_unknown), "")

  // Single-path identifier
  PeptideIdentification pep_a;
  pep_a.setIdentifier("ID_A");
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(pep_a), "/data/runA.mzML")
}
END_SECTION

START_SECTION([EXTRA] getPrimaryMSRunPath legacy base_name fallback)
{
  IdentifierMSRunMapper m(makeProtIds());

  // Backwards compatibility: a PeptideIdentification whose identifier is not in the mapping
  // but which carries the (deprecated) "base_name" meta value resolves to that base name.
  PeptideIdentification pep_legacy;
  pep_legacy.setIdentifier("DOES_NOT_EXIST");
  pep_legacy.setMetaValue(Constants::UserParam::BASE_NAME, "legacy_run.mzML");
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(pep_legacy), "legacy_run.mzML")

  // The primary MS run path from the mapping takes precedence over a legacy base_name.
  PeptideIdentification pep_both;
  pep_both.setIdentifier("ID_A");
  pep_both.setMetaValue(Constants::UserParam::BASE_NAME, "legacy_run.mzML");
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(pep_both), "/data/runA.mzML")
}
END_SECTION

START_SECTION([EXTRA] create() leaves a complete forward map even when it throws on duplicate paths)
{
  // Runs A and D have distinct paths; B and C collide on the same path. The collision
  // makes create() throw, but the forward map (identifier -> path) must still be fully
  // populated for *all* runs - including D, which is listed *after* the duplicate - so
  // that getPrimaryMSRunPath() does not silently drop run names in an order-dependent way
  // for callers (e.g. TextExporter's USI export) that catch the exception and continue.
  ProteinIdentification pA, pB, pC, pD;
  pA.setIdentifier("ID_A");
  pA.setPrimaryMSRunPath(StringList{"/data/runA.mzML"});
  pB.setIdentifier("ID_B");
  pB.setPrimaryMSRunPath(StringList{"/data/dup.mzML"});
  pC.setIdentifier("ID_C");
  pC.setPrimaryMSRunPath(StringList{"/data/dup.mzML"});
  pD.setIdentifier("ID_D");
  pD.setPrimaryMSRunPath(StringList{"/data/runD.mzML"});

  IdentifierMSRunMapper m;
  TEST_EXCEPTION(Exception::InvalidValue, m.create(std::vector<ProteinIdentification>{pA, pB, pC, pD}))

  // Despite the throw, every identifier still resolves through the forward map.
  PeptideIdentification pepA; pepA.setIdentifier("ID_A");
  PeptideIdentification pepD; pepD.setIdentifier("ID_D");
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(pepA), "/data/runA.mzML")
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(pepD), "/data/runD.mzML") // run listed after the duplicate
}
END_SECTION

START_SECTION((void validateMergeIndex(const PeptideIdentification& pep_id, size_t psm_index) const))
{
  // getPrimaryMSRunPath() falls back to the run's FIRST file when it cannot resolve an index,
  // which silently mislabels every PSM of a merged run. Callers that key on the result must be
  // able to refuse such input instead; this is that check.
  ProteinIdentification merged;
  merged.setIdentifier("MERGED");
  merged.setPrimaryMSRunPath(StringList{"/data/runA.mzML", "/data/runB.mzML"});
  ProteinIdentification single;
  single.setIdentifier("SINGLE");
  single.setPrimaryMSRunPath(StringList{"/data/only.mzML"});
  ProteinIdentification none;
  none.setIdentifier("NONE");

  IdentifierMSRunMapper m;
  m.create(std::vector<ProteinIdentification>{merged, single, none});

  // Accepted: a usable index into the merged run.
  PeptideIdentification ok;
  ok.setIdentifier("MERGED");
  ok.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, 1);
  m.validateMergeIndex(ok, 0);
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(ok), "/data/runB.mzML")

  // Exempt: unmerged runs have nothing to disambiguate, so a missing index is not an error.
  // Resolution is unambiguous either way, which is what the exemption rests on.
  PeptideIdentification one_path;
  one_path.setIdentifier("SINGLE");
  m.validateMergeIndex(one_path, 0);
  TEST_STRING_EQUAL(m.getPrimaryMSRunPath(one_path), "/data/only.mzML")

  PeptideIdentification no_path;
  no_path.setIdentifier("NONE");
  m.validateMergeIndex(no_path, 0);
  TEST_TRUE(m.getPrimaryMSRunPath(no_path).empty())

  // Refused: missing, wrongly typed, negative, and out-of-range indices.
  PeptideIdentification missing;
  missing.setIdentifier("MERGED");
  TEST_EXCEPTION(Exception::MissingInformation, m.validateMergeIndex(missing, 0))

  PeptideIdentification wrong_type;
  wrong_type.setIdentifier("MERGED");
  wrong_type.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, "1");
  TEST_EXCEPTION(Exception::MissingInformation, m.validateMergeIndex(wrong_type, 0))

  // Negative must surface as MissingInformation, not as the ConversionError that an
  // unsigned cast of the DataValue would raise.
  PeptideIdentification negative;
  negative.setIdentifier("MERGED");
  negative.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, -1);
  TEST_EXCEPTION(Exception::MissingInformation, m.validateMergeIndex(negative, 0))

  PeptideIdentification too_big;
  too_big.setIdentifier("MERGED");
  too_big.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, 2);
  TEST_EXCEPTION(Exception::MissingInformation, m.validateMergeIndex(too_big, 0))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
