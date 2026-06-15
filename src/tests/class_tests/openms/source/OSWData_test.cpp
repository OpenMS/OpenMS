// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/OSWData.h>
///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>

using namespace OpenMS;
using namespace std;

namespace
{
  // a small OSWData with two transitions and one protein referencing both
  OSWData makeData()
  {
    OSWData d;
    d.addTransition(OSWTransition("y5/-0.002", 100, 500.25f, 'y', false));
    d.addTransition(OSWTransition("b3", 101, 300.10f, 'b', true));

    std::vector<UInt32> tids{100, 101};
    OSWPeakGroup pg(120.5f, 118.0f, 123.0f, 1.5f, std::move(tids), 0.01f);
    std::vector<OSWPeakGroup> feats;
    feats.push_back(pg);
    OSWPeptidePrecursor pc("PEPTIDEK", 2, false, 450.7f, std::move(feats));
    std::vector<OSWPeptidePrecursor> peps;
    peps.push_back(pc);
    d.addProtein(OSWProtein("PROT1", 7, std::move(peps)));
    return d;
  }
}

START_TEST(OSWData, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] OSWHierarchy::LevelName)
{
  TEST_STRING_EQUAL(OSWHierarchy::LevelName[OSWHierarchy::PROTEIN], "protein")
  TEST_STRING_EQUAL(OSWHierarchy::LevelName[OSWHierarchy::PEPTIDE], "peptide")
  TEST_STRING_EQUAL(OSWHierarchy::LevelName[OSWHierarchy::FEATURE], "feature/peakgroup")
  TEST_STRING_EQUAL(OSWHierarchy::LevelName[OSWHierarchy::TRANSITION], "transition")
}
END_SECTION

START_SECTION([EXTRA] OSWIndexTrace::isSet)
{
  OSWIndexTrace t;
  TEST_EQUAL(t.isSet(), false) // default-constructed -> not pointing anywhere
  t.idx_prot = 0;
  t.lowest = OSWHierarchy::PROTEIN;
  TEST_EQUAL(t.isSet(), true)
}
END_SECTION

START_SECTION([EXTRA] OSWTransition value type)
{
  OSWTransition tr("y5/-0.002", 100, 500.25f, 'y', false);
  TEST_STRING_EQUAL(tr.getAnnotation(), "y5/-0.002")
  TEST_EQUAL(tr.getID(), 100)
  TEST_REAL_SIMILAR(tr.getProductMZ(), 500.25f)
  TEST_EQUAL(tr.getType(), 'y')
  TEST_EQUAL(tr.isDecoy(), false)
}
END_SECTION

START_SECTION([EXTRA] OSWPeakGroup value type)
{
  // default q-value is the "missing" sentinel
  OSWPeakGroup empty;
  TEST_REAL_SIMILAR(empty.getQValue(), OSWPeakGroup::QVALUE_MISSING)

  std::vector<UInt32> tids{100, 101};
  OSWPeakGroup pg(120.5f, 118.0f, 123.0f, 1.5f, std::move(tids), 0.01f);
  TEST_REAL_SIMILAR(pg.getRTExperimental(), 120.5f)
  TEST_REAL_SIMILAR(pg.getRTLeftWidth(), 118.0f)
  TEST_REAL_SIMILAR(pg.getRTRightWidth(), 123.0f)
  TEST_REAL_SIMILAR(pg.getRTDelta(), 1.5f)
  TEST_REAL_SIMILAR(pg.getQValue(), 0.01f)
  TEST_EQUAL(pg.getTransitionIDs().size(), 2)
  TEST_EQUAL(pg.getTransitionIDs()[0], 100)
  TEST_EQUAL(pg.getTransitionIDs()[1], 101)
}
END_SECTION

START_SECTION([EXTRA] OSWPeptidePrecursor value type)
{
  std::vector<OSWPeakGroup> feats;
  feats.push_back(OSWPeakGroup());
  OSWPeptidePrecursor pc("PEPTIDEK", 2, false, 450.7f, std::move(feats));
  TEST_STRING_EQUAL(pc.getSequence(), "PEPTIDEK")
  TEST_EQUAL(pc.getCharge(), 2)
  TEST_EQUAL(pc.isDecoy(), false)
  TEST_REAL_SIMILAR(pc.getPCMz(), 450.7f)
  TEST_EQUAL(pc.getFeatures().size(), 1)
}
END_SECTION

START_SECTION([EXTRA] OSWProtein value type)
{
  std::vector<OSWPeptidePrecursor> peps;
  peps.push_back(OSWPeptidePrecursor());
  OSWProtein prot("PROT1", 7, std::move(peps));
  TEST_STRING_EQUAL(prot.getAccession(), "PROT1")
  TEST_EQUAL(prot.getID(), 7)
  TEST_EQUAL(prot.getPeptidePrecursors().size(), 1)
}
END_SECTION

START_SECTION((void addTransition(const OSWTransition& tr)))
{
  OSWData d;
  OSWTransition tr("y5", 100, 500.0f, 'y', false);
  d.addTransition(tr); // const ref overload
  TEST_EQUAL(d.transitionCount(), 1)
  TEST_STRING_EQUAL(d.getTransition(100).getAnnotation(), "y5")
}
END_SECTION

START_SECTION((void addTransition(OSWTransition&& tr)))
{
  OSWData d;
  d.addTransition(OSWTransition("b3", 101, 300.0f, 'b', true)); // rvalue overload
  TEST_EQUAL(d.transitionCount(), 1)
  TEST_EQUAL(d.getTransition(101).isDecoy(), true)
}
END_SECTION

START_SECTION((Size transitionCount() const))
{
  TEST_EQUAL(makeData().transitionCount(), 2)
}
END_SECTION

START_SECTION((const OSWTransition& getTransition(const UInt32 id) const))
{
  OSWData d = makeData();
  TEST_STRING_EQUAL(d.getTransition(100).getAnnotation(), "y5/-0.002")
  TEST_EQUAL(d.getTransition(101).getType(), 'b')
}
END_SECTION

START_SECTION((const std::map<UInt32, OSWTransition>& getTransitions() const))
{
  OSWData d = makeData();
  TEST_EQUAL(d.getTransitions().size(), 2)
  TEST_EQUAL(d.getTransitions().count(100), 1)
  TEST_EQUAL(d.getTransitions().count(999), 0)
}
END_SECTION

START_SECTION((void addProtein(OSWProtein&& prot)))
{
  OSWData d = makeData();
  TEST_EQUAL(d.getProteins().size(), 1)

  // referential integrity: a protein referencing an unknown transition is rejected
  OSWData d2;
  d2.addTransition(OSWTransition("y5", 100, 500.0f, 'y', false));
  std::vector<UInt32> bad{999}; // not added
  OSWPeakGroup pg(1.0f, 0.0f, 2.0f, 0.0f, std::move(bad));
  std::vector<OSWPeakGroup> feats; feats.push_back(pg);
  OSWPeptidePrecursor pc("AAA", 1, false, 100.0f, std::move(feats));
  std::vector<OSWPeptidePrecursor> peps; peps.push_back(pc);
  TEST_EXCEPTION(Exception::Precondition, d2.addProtein(OSWProtein("P", 1, std::move(peps))))
}
END_SECTION

START_SECTION((const std::vector<OSWProtein>& getProteins() const))
{
  OSWData d = makeData();
  TEST_EQUAL(d.getProteins().size(), 1)
  TEST_STRING_EQUAL(d.getProteins()[0].getAccession(), "PROT1")
  // the full hierarchy is reachable through the accessors
  const auto& pcs = d.getProteins()[0].getPeptidePrecursors();
  TEST_EQUAL(pcs.size(), 1)
  TEST_STRING_EQUAL(pcs[0].getSequence(), "PEPTIDEK")
  TEST_EQUAL(pcs[0].getFeatures()[0].getTransitionIDs().size(), 2)
}
END_SECTION

START_SECTION((void setProtein(const Size index, OSWProtein&& protein)))
{
  OSWData d = makeData();
  std::vector<OSWPeptidePrecursor> peps; // stub protein, no transition refs
  d.setProtein(0, OSWProtein("REPLACED", 9, std::move(peps)));
  TEST_EQUAL(d.getProteins().size(), 1)
  TEST_STRING_EQUAL(d.getProteins()[0].getAccession(), "REPLACED")
}
END_SECTION

START_SECTION((void setSqlSourceFile(const std::string& filename)))
{
  OSWData d;
  d.setSqlSourceFile("/data/run.osw");
  TEST_STRING_EQUAL(d.getSqlSourceFile(), "/data/run.osw")
}
END_SECTION

START_SECTION((const std::string& getSqlSourceFile() const))
{
  OSWData d;
  TEST_STRING_EQUAL(d.getSqlSourceFile(), "")
  d.setSqlSourceFile("x.osw");
  TEST_STRING_EQUAL(d.getSqlSourceFile(), "x.osw")
}
END_SECTION

START_SECTION((void setRunID(const UInt64 run_id)))
{
  OSWData d;
  d.setRunID(123456789ULL);
  TEST_EQUAL(d.getRunID() == 123456789ULL, true)
}
END_SECTION

START_SECTION((UInt64 getRunID() const))
{
  OSWData d;
  d.setRunID(42);
  TEST_EQUAL(d.getRunID() == 42, true)
}
END_SECTION

START_SECTION((void clear()))
{
  OSWData d = makeData();
  d.clear();
  TEST_EQUAL(d.transitionCount(), 0)
  TEST_EQUAL(d.getProteins().size(), 0)
}
END_SECTION

START_SECTION((void clearProteins()))
{
  OSWData d = makeData();
  d.clearProteins();
  TEST_EQUAL(d.getProteins().size(), 0)
  TEST_EQUAL(d.transitionCount(), 2) // transitions are kept
}
END_SECTION

START_SECTION((void buildNativeIDResolver(const MSExperiment& chrom_traces)))
{
  OSWData d = makeData();
  d.setRunID(999); // will not match an unrelated/empty MSExperiment's run id
  MSExperiment exp;
  TEST_EXCEPTION(Exception::Precondition, d.buildNativeIDResolver(exp))
}
END_SECTION

START_SECTION((UInt fromNativeID(int transition_id) const))
{
  OSWData d = makeData();
  // without a prior buildNativeIDResolver() the mapping is empty -> throws
  TEST_EXCEPTION(Exception::InvalidValue, d.fromNativeID(100))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
