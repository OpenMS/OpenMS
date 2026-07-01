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

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/MRMFeature.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <algorithm>
#include <vector>

///////////////////////////

START_TEST(MRMFeatureAccessOpenMS, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// Helper: a Feature with RT, intensity and (optionally) a single convex hull
// whose two points are (95, 10) and (105, 20) -> X = RT, Y = intensity.
auto mkFeat = [](double rt, float inten, bool with_hull) -> Feature {
  Feature f;
  f.setRT(rt);
  f.setIntensity(inten);
  if (with_hull)
  {
    ConvexHull2D hull;
    ConvexHull2D::PointArrayType pts;
    pts.push_back(DPosition<2>(95.0, 10.0));
    pts.push_back(DPosition<2>(105.0, 20.0));
    hull.setHullPoints(pts);
    f.setConvexHulls(std::vector<ConvexHull2D> {hull});
  }
  return f;
};

// ---------------------------- FeatureOpenMS ----------------------------

FeatureOpenMS* ptr = nullptr;
FeatureOpenMS* null_ptr = nullptr;
Feature feat_for_ptr = mkFeat(100.0, 500.0f, true);

START_SECTION(([FeatureOpenMS] FeatureOpenMS(const Feature& feature)))
{
  ptr = new FeatureOpenMS(feat_for_ptr);
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(([FeatureOpenMS] ~FeatureOpenMS()))
{
  delete ptr;
}
END_SECTION

START_SECTION(([FeatureOpenMS] double getRT() const))
{
  Feature f = mkFeat(100.0, 500.0f, true);
  FeatureOpenMS fo(f);
  TEST_REAL_SIMILAR(fo.getRT(), 100.0)
}
END_SECTION

START_SECTION(([FeatureOpenMS] float getIntensity() const))
{
  Feature f = mkFeat(100.0, 500.0f, true);
  FeatureOpenMS fo(f);
  TEST_REAL_SIMILAR(fo.getIntensity(), 500.0)
}
END_SECTION

START_SECTION(([FeatureOpenMS] void getRT(std::vector<double>& rt) const))
{
  Feature f = mkFeat(100.0, 500.0f, true);
  FeatureOpenMS fo(f);
  std::vector<double> rt;
  fo.getRT(rt); // X coordinates of the single convex hull's points
  TEST_EQUAL(rt.size(), 2)
  TEST_REAL_SIMILAR(rt[0], 95.0)
  TEST_REAL_SIMILAR(rt[1], 105.0)

  // a feature without exactly one convex hull cannot provide the data points
  Feature no_hull = mkFeat(50.0, 5.0f, false);
  FeatureOpenMS fo_nohull(no_hull);
  std::vector<double> tmp;
  TEST_EXCEPTION(Exception::IllegalArgument, fo_nohull.getRT(tmp))
}
END_SECTION

START_SECTION(([FeatureOpenMS] void getIntensity(std::vector<double>& intens) const))
{
  Feature f = mkFeat(100.0, 500.0f, true);
  FeatureOpenMS fo(f);
  std::vector<double> intens;
  fo.getIntensity(intens); // Y coordinates of the single convex hull's points
  TEST_EQUAL(intens.size(), 2)
  TEST_REAL_SIMILAR(intens[0], 10.0)
  TEST_REAL_SIMILAR(intens[1], 20.0)

  Feature no_hull = mkFeat(50.0, 5.0f, false);
  FeatureOpenMS fo_nohull(no_hull);
  std::vector<double> tmp;
  TEST_EXCEPTION(Exception::IllegalArgument, fo_nohull.getIntensity(tmp))
}
END_SECTION

// ---------------------------- MRMFeatureOpenMS ----------------------------

START_SECTION(([MRMFeatureOpenMS] MRMFeatureOpenMS(MRMFeature & mrmfeature)))
{
  MRMFeature mrm;
  mrm.setRT(200.0);
  mrm.setIntensity(1000.0f);
  mrm.addFeature(mkFeat(198.0, 300.0f, true), "tr1");
  mrm.addFeature(mkFeat(202.0, 700.0f, true), "tr2");

  MRMFeatureOpenMS* mptr = new MRMFeatureOpenMS(mrm);
  TEST_NOT_EQUAL(mptr, null_ptr)
  delete mptr;
}
END_SECTION

START_SECTION(([MRMFeatureOpenMS] double getRT() const))
{
  MRMFeature mrm;
  mrm.setRT(200.0);
  mrm.addFeature(mkFeat(198.0, 300.0f, true), "tr1");
  MRMFeatureOpenMS mo(mrm);
  TEST_REAL_SIMILAR(mo.getRT(), 200.0)
}
END_SECTION

START_SECTION(([MRMFeatureOpenMS] float getIntensity() const))
{
  MRMFeature mrm;
  mrm.setIntensity(1000.0f);
  mrm.addFeature(mkFeat(198.0, 300.0f, true), "tr1");
  MRMFeatureOpenMS mo(mrm);
  TEST_REAL_SIMILAR(mo.getIntensity(), 1000.0)
}
END_SECTION

START_SECTION(([MRMFeatureOpenMS] size_t size() const))
{
  MRMFeature mrm;
  mrm.addFeature(mkFeat(198.0, 300.0f, true), "tr1");
  mrm.addFeature(mkFeat(202.0, 700.0f, true), "tr2");
  MRMFeatureOpenMS mo(mrm);
  TEST_EQUAL(mo.size(), 2)
}
END_SECTION

START_SECTION(([MRMFeatureOpenMS] double getMetaValue(std::string name) const))
{
  MRMFeature mrm;
  mrm.setMetaValue("score", 0.75);
  mrm.addFeature(mkFeat(198.0, 300.0f, true), "tr1");
  MRMFeatureOpenMS mo(mrm);
  TEST_REAL_SIMILAR(mo.getMetaValue("score"), 0.75)
}
END_SECTION

START_SECTION(([MRMFeatureOpenMS] std::vector<std::string> getNativeIDs() const))
{
  MRMFeature mrm;
  mrm.addFeature(mkFeat(198.0, 300.0f, true), "tr1");
  mrm.addFeature(mkFeat(202.0, 700.0f, true), "tr2");
  MRMFeatureOpenMS mo(mrm);

  std::vector<std::string> ids = mo.getNativeIDs();
  TEST_EQUAL(ids.size(), 2)
  TEST_EQUAL(std::find(ids.begin(), ids.end(), "tr1") != ids.end(), true)
  TEST_EQUAL(std::find(ids.begin(), ids.end(), "tr2") != ids.end(), true)
}
END_SECTION

START_SECTION(([MRMFeatureOpenMS] std::shared_ptr<OpenSwath::IFeature> getFeature(std::string nativeID)))
{
  MRMFeature mrm;
  mrm.addFeature(mkFeat(198.0, 300.0f, true), "tr1");
  mrm.addFeature(mkFeat(202.0, 700.0f, true), "tr2");
  MRMFeatureOpenMS mo(mrm);

  std::shared_ptr<OpenSwath::IFeature> f1 = mo.getFeature("tr1");
  TEST_REAL_SIMILAR(f1->getRT(), 198.0)
  TEST_REAL_SIMILAR(f1->getIntensity(), 300.0)

  std::shared_ptr<OpenSwath::IFeature> f2 = mo.getFeature("tr2");
  TEST_REAL_SIMILAR(f2->getIntensity(), 700.0)
}
END_SECTION

// ---------------------------- TransitionGroupOpenMS ----------------------------

typedef MRMTransitionGroup<MSChromatogram, ReactionMonitoringTransition> TransitionGroupType;

// Helper: build an MRMTransitionGroup with two transitions (library intensities)
// and two chromatograms (native ids).
auto mkTransitionGroup = []() -> TransitionGroupType {
  TransitionGroupType tg;
  ReactionMonitoringTransition t1;
  t1.setLibraryIntensity(0.8);
  t1.setNativeID("tr1");
  ReactionMonitoringTransition t2;
  t2.setLibraryIntensity(0.3);
  t2.setNativeID("tr2");
  tg.addTransition(t1, "tr1");
  tg.addTransition(t2, "tr2");
  MSChromatogram c1;
  c1.setNativeID("tr1");
  MSChromatogram c2;
  c2.setNativeID("tr2");
  tg.addChromatogram(c1, "tr1");
  tg.addChromatogram(c2, "tr2");
  return tg;
};

START_SECTION(([TransitionGroupOpenMS] std::size_t size() const))
{
  TransitionGroupType tg = mkTransitionGroup();
  TransitionGroupOpenMS<MSChromatogram, ReactionMonitoringTransition> tgo(tg);
  TEST_EQUAL(tgo.size(), 2) // number of chromatograms
}
END_SECTION

START_SECTION(([TransitionGroupOpenMS] std::vector<std::string> getNativeIDs() const))
{
  TransitionGroupType tg = mkTransitionGroup();
  TransitionGroupOpenMS<MSChromatogram, ReactionMonitoringTransition> tgo(tg);
  std::vector<std::string> ids = tgo.getNativeIDs(); // taken from the chromatogram native ids
  TEST_EQUAL(ids.size(), 2)
  TEST_STRING_EQUAL(ids[0], "tr1")
  TEST_STRING_EQUAL(ids[1], "tr2")
}
END_SECTION

START_SECTION(([TransitionGroupOpenMS] void getLibraryIntensities(std::vector<double>& intensities) const))
{
  TransitionGroupType tg = mkTransitionGroup();
  TransitionGroupOpenMS<MSChromatogram, ReactionMonitoringTransition> tgo(tg);
  std::vector<double> libs;
  tgo.getLibraryIntensities(libs);
  TEST_EQUAL(libs.size(), 2)
  TEST_REAL_SIMILAR(libs[0], 0.8)
  TEST_REAL_SIMILAR(libs[1], 0.3)
}
END_SECTION

// ---------------------------- SignalToNoiseOpenMS ----------------------------

START_SECTION(([SignalToNoiseOpenMS] double getValueAtRT(double RT)))
{
  // a flat baseline (intensity 10) with one sharp peak (intensity 1000) at RT 15
  MSChromatogram chrom;
  for (int i = 0; i < 30; ++i)
  {
    ChromatogramPeak p;
    p.setRT(i * 1.0);
    p.setIntensity(10.0);
    chrom.push_back(p);
  }
  chrom[15].setIntensity(1000.0);

  SignalToNoiseOpenMS<MSChromatogram> sn(chrom, 30.0, 30, false);
  double sn_apex = sn.getValueAtRT(15.0);
  double sn_base = sn.getValueAtRT(2.0);

  // the peak point has a far higher S/N than the baseline
  TEST_EQUAL(sn_apex > sn_base, true)
  // S/N is linear in intensity (shared noise estimate): apex (1000) is 100x baseline (10)
  TEST_REAL_SIMILAR(sn_apex / sn_base, 100.0)

  // an empty chromatogram yields the -1 sentinel
  MSChromatogram empty;
  SignalToNoiseOpenMS<MSChromatogram> sn_empty(empty, 30.0, 30, false);
  TEST_REAL_SIMILAR(sn_empty.getValueAtRT(5.0), -1.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
