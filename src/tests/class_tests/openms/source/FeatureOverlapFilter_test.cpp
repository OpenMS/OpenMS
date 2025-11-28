// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/PROCESSING/FEATURE/FeatureOverlapFilter.h>
#include <OpenMS/CONCEPT/Constants.h>

///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

// Helper function to create a simple feature with position, intensity, and charge
Feature createTestFeature(double rt, double mz, double intensity, int charge = 2)
{
  Feature f;
  f.setRT(rt);
  f.setMZ(mz);
  f.setIntensity(intensity);
  f.setCharge(charge);
  // Add a simple convex hull for the feature
  std::vector<ConvexHull2D> hulls(1);
  hulls[0].addPoint(DPosition<2>(rt - 1.0, mz - 0.01));
  hulls[0].addPoint(DPosition<2>(rt + 1.0, mz - 0.01));
  hulls[0].addPoint(DPosition<2>(rt + 1.0, mz + 0.01));
  hulls[0].addPoint(DPosition<2>(rt - 1.0, mz + 0.01));
  f.setConvexHulls(hulls);
  return f;
}

/////////////////////////////////////////////////////////////

START_TEST(FeatureOverlapFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION((Filter FeatureMap))
  //feature with convex hulls
  Feature feature1;
  feature1.getPosition()[0] = 5.25;
  feature1.getPosition()[1] = 1.5;
  feature1.setIntensity(0.5f);
  feature1.setOverallQuality(8);
  std::vector< ConvexHull2D > hulls(1);
  hulls[0].addPoint(DPosition<2>(-1.0,2.0));
  hulls[0].addPoint(DPosition<2>(4.0,1.2));
  hulls[0].addPoint(DPosition<2>(5.0,3.123));
  feature1.setConvexHulls(hulls);

  Feature feature2;
  feature2.getPosition()[0] = 5.25;
  feature2.getPosition()[1] = 1.5;
  feature2.setIntensity(0.5f);
  feature2.setOverallQuality(10);
  std::vector< ConvexHull2D > hulls2(1);
  hulls2[0].addPoint(DPosition<2>(-1.0,2.0));
  hulls2[0].addPoint(DPosition<2>(4.0,1.2));
  hulls2[0].addPoint(DPosition<2>(5.5,3.123));
  feature2.setConvexHulls(hulls2);

  Feature feature3;
  feature3.getPosition()[0] = 5.25;
  feature3.getPosition()[1] = 1.5;
  feature3.setIntensity(0.5f);
  feature3.setOverallQuality(7);
  std::vector< ConvexHull2D > hulls3(1);
  hulls3[0].addPoint(DPosition<2>(4.5,2.0));
  hulls3[0].addPoint(DPosition<2>(10,1.2));
  hulls3[0].addPoint(DPosition<2>(10,3.123));
  feature3.setConvexHulls(hulls3);

  Feature feature4;
  feature4.getPosition()[0] = 20.;
  feature4.getPosition()[1] = 10.;
  feature4.setIntensity(0.5f);
  feature4.setOverallQuality(7);
  std::vector< ConvexHull2D > hulls4(1);
  hulls4[0].addPoint(DPosition<2>(20,5));
  hulls4[0].addPoint(DPosition<2>(22,10));
  hulls4[0].addPoint(DPosition<2>(22,14));
  feature4.setConvexHulls(hulls4);

  Feature feature5;
  feature5.getPosition()[0] = 20.;
  feature5.getPosition()[1] = 11.;
  feature5.setIntensity(0.5f);
  feature5.setOverallQuality(0.);
  std::vector< ConvexHull2D > hulls5(1);
  hulls5[0].addPoint(DPosition<2>(20,12.));
  hulls5[0].addPoint(DPosition<2>(21,16.));
  hulls5[0].addPoint(DPosition<2>(21,18.));
  feature5.setConvexHulls(hulls5);

  FeatureMap fmap;
  fmap.emplace_back(feature1);
  fmap.emplace_back(feature2);
  fmap.emplace_back(feature3);
  fmap.emplace_back(feature4);
  fmap.emplace_back(feature5);

  fmap.updateRanges();
  for (auto& f : fmap)
  {
    f.ensureUniqueId();
  }
  
  FeatureOverlapFilter::filter(fmap, 
    [](const Feature& left, const Feature& right){ return left.getOverallQuality() > right.getOverallQuality(); }, 
    [](const Feature&, const Feature&) { return true; },
    false);

  TEST_EQUAL(fmap[0].getOverallQuality(), 10)
  TEST_EQUAL(fmap[1].getOverallQuality(), 7)

END_SECTION

START_SECTION(mergeOverlappingFeatures - basic merging with SUM intensity)
{
  FeatureMap fmap;

  // Two overlapping features (within 5.0 RT and 0.05 mz tolerance)
  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  Feature f2 = createTestFeature(102.0, 500.02, 500.0, 2);  // within tolerance
  // One non-overlapping feature
  Feature f3 = createTestFeature(200.0, 600.0, 800.0, 2);

  fmap.push_back(f1);
  fmap.push_back(f2);
  fmap.push_back(f3);

  for (auto& f : fmap) f.ensureUniqueId();

  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, false,
                                                  MergeIntensityMode::SUM, true);

  // Should have 2 features: one merged (f1+f2) and one separate (f3)
  TEST_EQUAL(fmap.size(), 2)

  // Find the merged feature (higher intensity, should be sum of 1000 + 500 = 1500)
  bool found_merged = false;
  bool found_separate = false;
  for (const auto& f : fmap)
  {
    if (std::abs(f.getIntensity() - 1500.0) < 0.01)
    {
      found_merged = true;
      // Check meta values were written
      TEST_EQUAL(f.metaValueExists("merged_centroid_rts"), true)
      TEST_EQUAL(f.metaValueExists("merged_centroid_mzs"), true)
      std::vector<double> merged_rts = f.getMetaValue("merged_centroid_rts");
      TEST_EQUAL(merged_rts.size(), 2)
    }
    if (std::abs(f.getIntensity() - 800.0) < 0.01)
    {
      found_separate = true;
      // Non-merged feature should not have merge meta values
      TEST_EQUAL(f.metaValueExists("merged_centroid_rts"), false)
    }
  }
  TEST_EQUAL(found_merged, true)
  TEST_EQUAL(found_separate, true)
}
END_SECTION

START_SECTION(mergeOverlappingFeatures - MAX intensity mode)
{
  FeatureMap fmap;

  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  Feature f2 = createTestFeature(102.0, 500.02, 500.0, 2);

  fmap.push_back(f1);
  fmap.push_back(f2);

  for (auto& f : fmap) f.ensureUniqueId();

  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, false,
                                                  MergeIntensityMode::MAX, true);

  TEST_EQUAL(fmap.size(), 1)
  // MAX mode should keep 1000.0 (the maximum)
  TEST_REAL_SIMILAR(fmap[0].getIntensity(), 1000.0)
}
END_SECTION

START_SECTION(mergeOverlappingFeatures - require_same_charge)
{
  FeatureMap fmap;

  // Two features at same position but different charges
  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  Feature f2 = createTestFeature(101.0, 500.01, 500.0, 3);  // different charge

  fmap.push_back(f1);
  fmap.push_back(f2);

  for (auto& f : fmap) f.ensureUniqueId();

  // With require_same_charge=true, features should NOT merge
  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, false,
                                                  MergeIntensityMode::SUM, true);

  TEST_EQUAL(fmap.size(), 2)

  // Reset and test with require_same_charge=false
  fmap.clear();
  f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  f2 = createTestFeature(101.0, 500.01, 500.0, 3);
  fmap.push_back(f1);
  fmap.push_back(f2);
  for (auto& f : fmap) f.ensureUniqueId();

  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, false, false,
                                                  MergeIntensityMode::SUM, true);

  // With require_same_charge=false, features SHOULD merge
  TEST_EQUAL(fmap.size(), 1)
  TEST_REAL_SIMILAR(fmap[0].getIntensity(), 1500.0)
}
END_SECTION

START_SECTION(mergeOverlappingFeatures - require_same_im with FAIMS_CV)
{
  FeatureMap fmap;

  // Two features with different FAIMS CV values
  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  f1.setMetaValue(Constants::UserParam::FAIMS_CV, -45.0);
  Feature f2 = createTestFeature(101.0, 500.01, 500.0, 2);
  f2.setMetaValue(Constants::UserParam::FAIMS_CV, -60.0);  // different CV

  fmap.push_back(f1);
  fmap.push_back(f2);

  for (auto& f : fmap) f.ensureUniqueId();

  // With require_same_im=true, features should NOT merge (different CVs)
  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, true,
                                                  MergeIntensityMode::SUM, true);

  TEST_EQUAL(fmap.size(), 2)

  // Reset and test with require_same_im=false
  fmap.clear();
  f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  f1.setMetaValue(Constants::UserParam::FAIMS_CV, -45.0);
  f2 = createTestFeature(101.0, 500.01, 500.0, 2);
  f2.setMetaValue(Constants::UserParam::FAIMS_CV, -60.0);
  fmap.push_back(f1);
  fmap.push_back(f2);
  for (auto& f : fmap) f.ensureUniqueId();

  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, false,
                                                  MergeIntensityMode::SUM, true);

  // With require_same_im=false, features SHOULD merge
  TEST_EQUAL(fmap.size(), 1)
  TEST_REAL_SIMILAR(fmap[0].getIntensity(), 1500.0)

  // Check FAIMS CV values were collected
  TEST_EQUAL(fmap[0].metaValueExists("merged_centroid_IMs"), true)
  std::vector<double> merged_ims = fmap[0].getMetaValue("merged_centroid_IMs");
  TEST_EQUAL(merged_ims.size(), 2)
  TEST_EQUAL(fmap[0].getMetaValue("FAIMS_merge_count"), 2)
}
END_SECTION

START_SECTION(mergeOverlappingFeatures - require_same_im with same FAIMS_CV)
{
  FeatureMap fmap;

  // Two features with same FAIMS CV values
  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  f1.setMetaValue(Constants::UserParam::FAIMS_CV, -45.0);
  Feature f2 = createTestFeature(101.0, 500.01, 500.0, 2);
  f2.setMetaValue(Constants::UserParam::FAIMS_CV, -45.0);  // same CV

  fmap.push_back(f1);
  fmap.push_back(f2);

  for (auto& f : fmap) f.ensureUniqueId();

  // With require_same_im=true, features SHOULD merge (same CVs)
  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, true,
                                                  MergeIntensityMode::SUM, true);

  TEST_EQUAL(fmap.size(), 1)
  TEST_REAL_SIMILAR(fmap[0].getIntensity(), 1500.0)
}
END_SECTION

START_SECTION(mergeOverlappingFeatures - features without FAIMS_CV)
{
  FeatureMap fmap;

  // Two features without FAIMS CV (should merge when require_same_im=true)
  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  Feature f2 = createTestFeature(101.0, 500.01, 500.0, 2);

  fmap.push_back(f1);
  fmap.push_back(f2);

  for (auto& f : fmap) f.ensureUniqueId();

  // Both without FAIMS_CV should be treated as same group
  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, true,
                                                  MergeIntensityMode::SUM, true);

  TEST_EQUAL(fmap.size(), 1)
  TEST_REAL_SIMILAR(fmap[0].getIntensity(), 1500.0)
  // No FAIMS CV meta values should be present
  TEST_EQUAL(fmap[0].metaValueExists("merged_centroid_IMs"), false)
  TEST_EQUAL(fmap[0].metaValueExists("FAIMS_merge_count"), false)
}
END_SECTION

START_SECTION(mergeOverlappingFeatures - mixed FAIMS_CV presence with require_same_im)
{
  FeatureMap fmap;

  // One feature with FAIMS CV, one without
  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  f1.setMetaValue(Constants::UserParam::FAIMS_CV, -45.0);
  Feature f2 = createTestFeature(101.0, 500.01, 500.0, 2);
  // f2 has no FAIMS_CV

  fmap.push_back(f1);
  fmap.push_back(f2);

  for (auto& f : fmap) f.ensureUniqueId();

  // With require_same_im=true, features should NOT merge (one has CV, one doesn't)
  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, true,
                                                  MergeIntensityMode::SUM, true);

  TEST_EQUAL(fmap.size(), 2)
}
END_SECTION

START_SECTION(mergeOverlappingFeatures - write_meta_values=false)
{
  FeatureMap fmap;

  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  f1.setMetaValue(Constants::UserParam::FAIMS_CV, -45.0);
  Feature f2 = createTestFeature(101.0, 500.01, 500.0, 2);
  f2.setMetaValue(Constants::UserParam::FAIMS_CV, -60.0);

  fmap.push_back(f1);
  fmap.push_back(f2);

  for (auto& f : fmap) f.ensureUniqueId();

  // With write_meta_values=false, no merge tracking meta values should be written
  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, false,
                                                  MergeIntensityMode::SUM, false);

  TEST_EQUAL(fmap.size(), 1)
  TEST_REAL_SIMILAR(fmap[0].getIntensity(), 1500.0)
  // No meta values should be present
  TEST_EQUAL(fmap[0].metaValueExists("merged_centroid_rts"), false)
  TEST_EQUAL(fmap[0].metaValueExists("merged_centroid_mzs"), false)
  TEST_EQUAL(fmap[0].metaValueExists("merged_centroid_IMs"), false)
  TEST_EQUAL(fmap[0].metaValueExists("FAIMS_merge_count"), false)
}
END_SECTION

START_SECTION(mergeOverlappingFeatures - no merge when outside tolerance)
{
  FeatureMap fmap;

  // Two features outside the tolerance
  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  Feature f2 = createTestFeature(110.0, 500.0, 500.0, 2);  // 10 seconds apart (> 5.0 tolerance)

  fmap.push_back(f1);
  fmap.push_back(f2);

  for (auto& f : fmap) f.ensureUniqueId();

  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, false,
                                                  MergeIntensityMode::SUM, true);

  // Features should NOT merge (outside RT tolerance)
  TEST_EQUAL(fmap.size(), 2)
}
END_SECTION

START_SECTION(mergeOverlappingFeatures - multiple features merging)
{
  FeatureMap fmap;

  // Three features that should all merge together
  Feature f1 = createTestFeature(100.0, 500.0, 1000.0, 2);
  Feature f2 = createTestFeature(101.0, 500.01, 500.0, 2);
  Feature f3 = createTestFeature(102.0, 500.02, 300.0, 2);

  fmap.push_back(f1);
  fmap.push_back(f2);
  fmap.push_back(f3);

  for (auto& f : fmap) f.ensureUniqueId();

  FeatureOverlapFilter::mergeOverlappingFeatures(fmap, 5.0, 0.05, true, false,
                                                  MergeIntensityMode::SUM, true);

  TEST_EQUAL(fmap.size(), 1)
  // All three intensities summed
  TEST_REAL_SIMILAR(fmap[0].getIntensity(), 1800.0)

  // Check all three positions were collected
  std::vector<double> merged_rts = fmap[0].getMetaValue("merged_centroid_rts");
  TEST_EQUAL(merged_rts.size(), 3)
}
END_SECTION

END_TEST
