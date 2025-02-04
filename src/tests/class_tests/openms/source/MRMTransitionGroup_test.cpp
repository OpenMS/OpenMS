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

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

typedef OpenMS::ReactionMonitoringTransition TransitionType;
typedef MRMTransitionGroup<MSChromatogram, TransitionType> MRMTransitionGroupType;

///////////////////////////

START_TEST(MRMTransitionGroup, "$Id$")

/////////////////////////////////////////////////////////////

MRMTransitionGroupType* ptr = nullptr;
MRMTransitionGroupType* nullPointer = nullptr;

START_SECTION(MRMTransitionGroup())
{
	ptr = new MRMTransitionGroupType();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMTransitionGroup())
{
  delete ptr;
}
END_SECTION


MSChromatogram chrom1;
MSChromatogram chrom2;
TransitionType trans1;
TransitionType trans2;
MRMFeature feature1;
MRMFeature feature2;

START_SECTION(MRMTransitionGroup(const MRMTransitionGroup &rhs))
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");

  MRMTransitionGroupType tmp(mrmtrgroup);
  TEST_EQUAL(mrmtrgroup.size(), tmp.size() )
}
END_SECTION

    
START_SECTION( MRMTransitionGroup& operator=(const MRMTransitionGroup &rhs) )
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");

  MRMTransitionGroupType tmp = mrmtrgroup;
  TEST_EQUAL(mrmtrgroup.size(), tmp.size() )
}
END_SECTION

START_SECTION (Size size() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.size(), 1)
  mrmtrgroup.addChromatogram(chrom2, "dummy2");
  TEST_EQUAL(mrmtrgroup.size(), 2)
}
END_SECTION

START_SECTION (  const String & getTransitionGroupID() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.setTransitionGroupID("some_id");
  TEST_EQUAL(mrmtrgroup.getTransitionGroupID(), "some_id")
}
END_SECTION

START_SECTION (  void setTransitionGroupID(const String & tr_gr_id))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION ( std::vector<TransitionType>& getTransitionsMuteable())
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addTransition(trans1, "dummy1");
  mrmtrgroup.addTransition(trans2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getTransitionsMuteable().size(), 2)
}
END_SECTION

START_SECTION ( void addTransition(const TransitionType &transition, String key))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION ( const TransitionType& getTransition(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  trans1.setLibraryIntensity(42);
  mrmtrgroup.addTransition(trans1, "dummy1");
  TEST_EQUAL(mrmtrgroup.getTransition("dummy1").getLibraryIntensity(), 42)
}
END_SECTION

START_SECTION (  const std::vector<TransitionType>& getTransitions() const )
{
  MRMTransitionGroupType mrmtrgroup;
  trans1.setLibraryIntensity(42);
  mrmtrgroup.addTransition(trans1, "dummy1");
  trans2.setLibraryIntensity(-2);
  mrmtrgroup.addTransition(trans2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getTransitions()[0].getLibraryIntensity(), 42)
  TEST_EQUAL(mrmtrgroup.getTransitions()[1].getLibraryIntensity(), -2)
}
END_SECTION

START_SECTION (  bool hasTransition(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addTransition(trans1, "dummy1");
  TEST_EQUAL(mrmtrgroup.hasTransition("dummy1"), true)
  TEST_EQUAL(mrmtrgroup.hasTransition("dummy2"), false)
}
END_SECTION

START_SECTION (  const std::vector<SpectrumType>& getChromatograms() const ) 
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getChromatograms().size(), 2)
}
END_SECTION

START_SECTION (  std::vector<SpectrumType>& getChromatograms())
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getChromatograms().size(), 2)
}
END_SECTION

START_SECTION (  void addChromatogram(SpectrumType &chromatogram, String key)) 
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION (  SpectrumType& getChromatogram(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  chrom1.setMetaValue("some_value", 1);
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.getChromatogram("dummy1").getMetaValue("some_value"), 1)
}
END_SECTION

START_SECTION (  bool hasChromatogram(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.hasChromatogram("dummy1"), true)
  TEST_EQUAL(mrmtrgroup.hasChromatogram("dummy2"), false)
}
END_SECTION

START_SECTION (  void addPrecusorChromatogram(SpectrumType &chromatogram, String key)) 
{
  // tested below
  NOT_TESTABLE
}
END_SECTION

START_SECTION (  SpectrumType& getPrecursorChromatogram(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  chrom1.setMetaValue("some_value", 1);
  mrmtrgroup.addPrecursorChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.getPrecursorChromatogram("dummy1").getMetaValue("some_value"), 1)

  // Add a few feature chromatograms and then add a precursor chromatogram -> it should still work
  mrmtrgroup.addChromatogram(chrom1, "feature1");
  mrmtrgroup.addChromatogram(chrom1, "feature2");
  mrmtrgroup.addChromatogram(chrom1, "feature3");
  mrmtrgroup.addPrecursorChromatogram(chrom1, "dummy2");
  TEST_EQUAL(mrmtrgroup.getPrecursorChromatogram("dummy2").getMetaValue("some_value"), 1)
}
END_SECTION

START_SECTION (  bool hasPrecursorChromatogram(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addPrecursorChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.hasPrecursorChromatogram("dummy1"), true)
  TEST_EQUAL(mrmtrgroup.hasPrecursorChromatogram("dummy2"), false)
}
END_SECTION

START_SECTION (  const std::vector<MRMFeature> & getFeatures() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addFeature(feature1);
  mrmtrgroup.addFeature(feature2);
  TEST_EQUAL(mrmtrgroup.getFeatures().size(), 2)
}
END_SECTION

START_SECTION (  std::vector<MRMFeature> & getFeaturesMuteable())
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addFeature(feature1);
  mrmtrgroup.addFeature(feature2);
  TEST_EQUAL(mrmtrgroup.getFeaturesMuteable().size(), 2)
}
END_SECTION

START_SECTION (  void addFeature(MRMFeature & feature))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION ( void getLibraryIntensity(std::vector<double> & result) const)
{
  TransitionType new_trans1;
  TransitionType new_trans2;
  MRMTransitionGroupType mrmtrgroup;
  new_trans1.setLibraryIntensity(3);
  new_trans2.setLibraryIntensity(-2);
  mrmtrgroup.addTransition(new_trans1, "dummy1");
  mrmtrgroup.addTransition(new_trans2, "dummy2");
  std::vector< double > result;
  mrmtrgroup.getLibraryIntensity(result);
  TEST_EQUAL(result.size(), 2)
  TEST_REAL_SIMILAR(result[0], 3)
  TEST_REAL_SIMILAR(result[1], 0)
}
END_SECTION

START_SECTION ( MRMTransitionGroup subset(std::vector<std::string> tr_ids))
{
  TransitionType new_trans1;
  TransitionType new_trans2;
  MRMTransitionGroupType mrmtrgroup, mrmtrgroupsub;
  new_trans1.setLibraryIntensity(3);
  new_trans1.setNativeID("new_trans1");
  new_trans1.setMetaValue("detecting_transition","true");
  new_trans2.setLibraryIntensity(-2);
  new_trans2.setNativeID("new_trans2");
  new_trans2.setMetaValue("detecting_transition","false");
  mrmtrgroup.addTransition(new_trans1, "new_trans1");
  mrmtrgroup.addTransition(new_trans2, "new_trans2");
  std::vector< std::string > transition_ids;
  transition_ids.push_back("new_trans1");

  std::vector< double > result;
  mrmtrgroupsub = mrmtrgroup.subset(transition_ids);
  mrmtrgroupsub.getLibraryIntensity(result);
  TEST_EQUAL(result.size(), 1)
  TEST_REAL_SIMILAR(result[0], 3)
}
END_SECTION

START_SECTION ( inline bool isInternallyConsistent() const)
{
  MRMTransitionGroupType mrmtrgroup;
  TEST_EQUAL(mrmtrgroup.isInternallyConsistent(), true)
}
END_SECTION

START_SECTION (inline bool chromatogramIdsMatch() const)
{
  
  {
    MRMTransitionGroupType mrmtrgroup;
    Chromatogram c;
    c.setNativeID("test");
    mrmtrgroup.addChromatogram(c, "test");

    TEST_EQUAL(mrmtrgroup.chromatogramIdsMatch(), true)
    mrmtrgroup.addChromatogram(c, "test2");
    TEST_EQUAL(mrmtrgroup.chromatogramIdsMatch(), false)
  }
  

  {
    MRMTransitionGroupType mrmtrgroup;
    Chromatogram c;
    c.setNativeID("test");
    mrmtrgroup.addPrecursorChromatogram(c, "test");

    TEST_EQUAL(mrmtrgroup.chromatogramIdsMatch(), true)
    mrmtrgroup.addPrecursorChromatogram(c, "test2");
    TEST_EQUAL(mrmtrgroup.chromatogramIdsMatch(), false)
  }

}
END_SECTION


START_SECTION ( MRMTransitionGroup subsetDependent(std::vector<std::string> tr_ids))
{
  TransitionType new_trans1;
  TransitionType new_trans2;
  MRMTransitionGroupType mrmtrgroup, mrmtrgroupsub;
  new_trans1.setLibraryIntensity(3);
  new_trans1.setNativeID("new_trans1");
  new_trans1.setMetaValue("detecting_transition","true");
  new_trans2.setLibraryIntensity(-2);
  new_trans2.setNativeID("new_trans2");
  new_trans2.setMetaValue("detecting_transition","false");
  mrmtrgroup.addTransition(new_trans1, "new_trans1");
  mrmtrgroup.addTransition(new_trans2, "new_trans2");
  std::vector< std::string > transition_ids;
  transition_ids.push_back("new_trans1");
  transition_ids.push_back("new_trans2");

  std::vector< double > result;
  mrmtrgroupsub = mrmtrgroup.subset(transition_ids);
  mrmtrgroupsub.getLibraryIntensity(result);
  TEST_EQUAL(result.size(), 2)
  TEST_REAL_SIMILAR(result[0], 3)
  TEST_REAL_SIMILAR(result[1], 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST

