// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright The OpenMS team, Eberhard Karls University Tübingen,
//  ETH Zürich and FU Berlin 2001-2012.
//  This software is released under a BSD license. For a full list of
//  authors, refer to the file AUTHORS. For full licensing conditions
//  refer to the file LICENSE.
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>

using namespace OpenMS;
using namespace std;

typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;
typedef MRMTransitionGroup<MSSpectrum, ChromatogramPeak, OpenMS::ReactionMonitoringTransition> MRMTransitionGroupType;
typedef OpenMS::ReactionMonitoringTransition TransitionType;

///////////////////////////

START_TEST(MRMTransitionGroup, "$Id$")

/////////////////////////////////////////////////////////////

MRMTransitionGroupType* ptr = 0;
MRMTransitionGroupType* nullPointer = 0;

START_SECTION(TransitionGroup())
{
	ptr = new MRMTransitionGroupType();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TransitionGroup())
{
  delete ptr;
}
END_SECTION

RichPeakChromatogram chrom1;
RichPeakChromatogram chrom2;
TransitionType trans1;
TransitionType trans2;
MRMFeature feature1;
MRMFeature feature2;

START_SECTION ( inline int size() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.size(), 1)
  //mrmtrgroup.addChromatogram(chrom2, "dummy1");
  //TEST_EQUAL(mrmtrgroup.size(), 1)
  mrmtrgroup.addChromatogram(chrom2, "dummy2");
  TEST_EQUAL(mrmtrgroup.size(), 2)
}
END_SECTION

START_SECTION ( inline const String & getTransitionGroupID() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.setTransitionGroupID("some_id");
  TEST_EQUAL(mrmtrgroup.getTransitionGroupID(), "some_id")
}
END_SECTION

START_SECTION ( inline void setTransitionGroupID(const String & tr_gr_id))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION ( inline const std::vector<ReactionMonitoringTransition> & getTransitions() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addTransition(trans1, "dummy1");
  mrmtrgroup.addTransition(trans2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getTransitions().size(), 2)
}
END_SECTION

START_SECTION ( inline std::vector<ReactionMonitoringTransition> & getTransitionsMuteable())
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addTransition(trans1, "dummy1");
  mrmtrgroup.addTransition(trans2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getTransitionsMuteable().size(), 2)
}
END_SECTION

START_SECTION ( inline void addTransition(const ReactionMonitoringTransition & transition, String key))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION ( inline const ReactionMonitoringTransition & getTransition(String key) )
{
  MRMTransitionGroupType mrmtrgroup;
  trans1.setLibraryIntensity(42);
  mrmtrgroup.addTransition(trans1, "dummy1");
  TEST_EQUAL(mrmtrgroup.getTransition("dummy1").getLibraryIntensity(), 42)
}
END_SECTION

START_SECTION ( inline bool hasTransition(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addTransition(trans1, "dummy1");
  TEST_EQUAL(mrmtrgroup.hasTransition("dummy1"), true)
  TEST_EQUAL(mrmtrgroup.hasTransition("dummy2"), false)
}
END_SECTION

START_SECTION ( inline const std::vector<SpectrumPeakType> & getChromatograms() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getChromatograms().size(), 2)
}
END_SECTION

START_SECTION ( inline std::vector<SpectrumPeakType> & getChromatograms())
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  mrmtrgroup.addChromatogram(chrom2, "dummy2");
  TEST_EQUAL(mrmtrgroup.getChromatograms().size(), 2)
}
END_SECTION

START_SECTION ( inline void addChromatogram(SpectrumPeakType & chromatogram, String key))
{
  // tested above
  NOT_TESTABLE
}
END_SECTION

START_SECTION ( inline SpectrumPeakType & getChromatogram(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  chrom1.setMetaValue("some_value", 1);
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.getChromatogram("dummy1").getMetaValue("some_value"), 1)
}
END_SECTION

START_SECTION ( inline bool hasChromatogram(String key))
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addChromatogram(chrom1, "dummy1");
  TEST_EQUAL(mrmtrgroup.hasChromatogram("dummy1"), true)
  TEST_EQUAL(mrmtrgroup.hasChromatogram("dummy2"), false)
}
END_SECTION

START_SECTION ( inline const std::vector<MRMFeature> & getFeatures() const)
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addFeature(feature1);
  mrmtrgroup.addFeature(feature2);
  TEST_EQUAL(mrmtrgroup.getFeatures().size(), 2)
}
END_SECTION

START_SECTION ( inline std::vector<MRMFeature> & getFeaturesMuteable())
{
  MRMTransitionGroupType mrmtrgroup;
  mrmtrgroup.addFeature(feature1);
  mrmtrgroup.addFeature(feature2);
  TEST_EQUAL(mrmtrgroup.getFeaturesMuteable().size(), 2)
}
END_SECTION

START_SECTION ( inline void addFeature(MRMFeature & feature))
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

/////////////////////////////////////////////////////////////
END_TEST
