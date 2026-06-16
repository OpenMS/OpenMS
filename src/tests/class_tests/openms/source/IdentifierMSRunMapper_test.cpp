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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
