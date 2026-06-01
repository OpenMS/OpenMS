// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: George Rosenberger, Hannes Roest, Witold Wolski $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <boost/assign/std/vector.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/MRMFeatureAccessOpenMS.h>
///////////////////////////

#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

namespace
{
  Feature makeFeature_(double rt, double intensity, const std::vector<std::pair<double, double>>& hull_points)
  {
    Feature feature;
    feature.setRT(rt);
    feature.setIntensity(intensity);

    ConvexHull2D hull;
    ConvexHull2D::PointArrayType points;
    for (const auto& hull_point : hull_points)
    {
      DPosition<2> position;
      position[0] = hull_point.first;
      position[1] = hull_point.second;
      points.push_back(position);
    }
    hull.setHullPoints(points);
    feature.getConvexHulls().push_back(hull);
    return feature;
  }

  MRMFeature makeMRMFeatureFixture_()
  {
    MRMFeature feature;
    feature.addFeature(makeFeature_(10.0, 100.0, {{9.9, 11.0}, {10.1, 12.0}}), "tr2");
    feature.addFeature(makeFeature_(11.0, 200.0, {{10.9, 21.0}, {11.1, 22.0}}), "tr1");
    feature.addPrecursorFeature(makeFeature_(12.0, 300.0, {{11.9, 31.0}, {12.1, 32.0}}), "pr1");
    feature.addPrecursorFeature(makeFeature_(13.0, 400.0, {{12.9, 41.0}, {13.1, 42.0}}), "pr0");
    return feature;
  }
}

START_TEST(MRMFeatureAccessOpenMS, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

//FeatureOpenMS
{
FeatureOpenMS* ptr = nullptr;
FeatureOpenMS* nullPointer = nullptr;

START_SECTION(FeatureOpenMS())
{
  Feature f;
  ptr = new FeatureOpenMS(f);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~FeatureOpenMS())
{
  delete ptr;
}
END_SECTION
}

//MRMFeatureOpenMS
{
MRMFeatureOpenMS* ptr = nullptr;
MRMFeatureOpenMS* nullPointer = nullptr;

START_SECTION(MRMFeatureOpenMS())
{
  MRMFeature f;
  ptr = new MRMFeatureOpenMS(f);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMFeatureOpenMS())
{
  delete ptr;
}
END_SECTION

START_SECTION([EXTRA] default lookup mode preserves native-id access and falls back from wrong index hints)
{
  MRMFeature feature = makeMRMFeatureFixture_();
  MRMFeatureOpenMS access(feature);

  const std::vector<std::string> native_ids = access.getNativeIDs();
  TEST_EQUAL(native_ids.size(), 2)
  TEST_EQUAL(native_ids[0], "tr1")
  TEST_EQUAL(native_ids[1], "tr2")

  const std::vector<std::string> precursor_ids = access.getPrecursorIDs();
  TEST_EQUAL(precursor_ids.size(), 2)
  TEST_EQUAL(precursor_ids[0], "pr0")
  TEST_EQUAL(precursor_ids[1], "pr1")

  TEST_REAL_SIMILAR(access.getFeatureIntensity("tr1", 0), 200.0)
  TEST_REAL_SIMILAR(access.getFeatureIntensity("tr2", 1), 100.0)
  TEST_REAL_SIMILAR(access.getFeatureIntensity("tr2", 0), 100.0)
}
END_SECTION

START_SECTION([EXTRA] cached-id constructor keeps requested id order for feature and precursor extraction)
{
  MRMFeature feature = makeMRMFeatureFixture_();
  std::vector<std::string> feature_ids{"tr2", "tr1"};
  std::vector<std::string> precursor_ids{"pr1", "pr0"};
  MRMFeatureOpenMS access(feature, feature_ids, precursor_ids);

  const std::vector<std::string> native_ids = access.getNativeIDs();
  TEST_EQUAL(native_ids.size(), 2)
  TEST_EQUAL(native_ids[0], "tr2")
  TEST_EQUAL(native_ids[1], "tr1")

  const std::vector<std::string> native_precursor_ids = access.getPrecursorIDs();
  TEST_EQUAL(native_precursor_ids.size(), 2)
  TEST_EQUAL(native_precursor_ids[0], "pr1")
  TEST_EQUAL(native_precursor_ids[1], "pr0")
  TEST_REAL_SIMILAR(access.getFeatureIntensity("tr2", 0), 100.0)
  TEST_REAL_SIMILAR(access.getFeatureIntensity("tr1", 1), 200.0)

  std::vector<std::vector<double>> intensities;
  access.getFeatureIntensities(feature_ids, intensities);
  TEST_EQUAL(intensities.size(), 2)
  TEST_EQUAL(intensities[0].size(), 2)
  TEST_REAL_SIMILAR(intensities[0][0], 11.0)
  TEST_REAL_SIMILAR(intensities[0][1], 12.0)
  TEST_REAL_SIMILAR(intensities[1][0], 21.0)
  TEST_REAL_SIMILAR(intensities[1][1], 22.0)

  std::vector<std::vector<double>> precursor_intensities;
  access.getPrecursorFeatureIntensities(precursor_ids, precursor_intensities);
  TEST_EQUAL(precursor_intensities.size(), 2)
  TEST_EQUAL(precursor_intensities[0].size(), 2)
  TEST_REAL_SIMILAR(precursor_intensities[0][0], 31.0)
  TEST_REAL_SIMILAR(precursor_intensities[0][1], 32.0)
  TEST_REAL_SIMILAR(precursor_intensities[1][0], 41.0)
  TEST_REAL_SIMILAR(precursor_intensities[1][1], 42.0)
}
END_SECTION

START_SECTION([EXTRA] aligned-order constructor uses requested aligned ids and tolerates wrong index hints)
{
  MRMFeature feature = makeMRMFeatureFixture_();
  std::vector<std::string> feature_ids{"tr2", "tr1"};
  std::vector<std::string> precursor_ids{"pr1", "pr0"};
  std::vector<String> feature_lookup_ids{"tr1", "tr2"};
  std::vector<String> precursor_lookup_ids{"pr1", "pr0"};
  MRMFeatureOpenMS access(feature, feature_ids, precursor_ids, feature_lookup_ids, precursor_lookup_ids);

  const std::vector<std::string> native_ids = access.getNativeIDs();
  TEST_EQUAL(native_ids.size(), 2)
  TEST_EQUAL(native_ids[0], "tr2")
  TEST_EQUAL(native_ids[1], "tr1")

  const std::vector<std::string> native_precursor_ids = access.getPrecursorIDs();
  TEST_EQUAL(native_precursor_ids.size(), 2)
  TEST_EQUAL(native_precursor_ids[0], "pr1")
  TEST_EQUAL(native_precursor_ids[1], "pr0")
  TEST_REAL_SIMILAR(access.getFeatureIntensity("tr2", 0), 100.0)
  TEST_REAL_SIMILAR(access.getFeatureIntensity("tr1", 1), 200.0)
  TEST_REAL_SIMILAR(access.getFeatureIntensity("tr2", 1), 100.0)

  std::shared_ptr<OpenSwath::IFeature> wrapped_feature = access.getFeature("tr1");
  TEST_REAL_SIMILAR(wrapped_feature->getIntensity(), 200.0)

  std::vector<std::vector<double>> intensities;
  access.getFeatureIntensities(feature_ids, intensities);
  TEST_EQUAL(intensities.size(), 2)
  TEST_REAL_SIMILAR(intensities[0][0], 11.0)
  TEST_REAL_SIMILAR(intensities[1][0], 21.0)
}
END_SECTION
}

//TransitionGroupOpenMS
{
TransitionGroupOpenMS <MSChromatogram, ReactionMonitoringTransition>* ptr = nullptr;
TransitionGroupOpenMS <MSChromatogram, ReactionMonitoringTransition>* nullPointer = nullptr;

START_SECTION(TransitionGroupOpenMS())
{
  MRMTransitionGroup <MSChromatogram, ReactionMonitoringTransition> trgroup;
  ptr = new TransitionGroupOpenMS < MSChromatogram, ReactionMonitoringTransition> (trgroup);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TransitionGroupOpenMS())
{
  delete ptr;
}
END_SECTION
}

//SignalToNoiseOpenMS
{
SignalToNoiseOpenMS<MSSpectrum>* ptr = nullptr;
SignalToNoiseOpenMS<MSSpectrum>* nullPointer = nullptr;

START_SECTION(SignalToNoiseOpenMS())
{
  OpenMS::MSSpectrum chromat;
  ptr = new SignalToNoiseOpenMS<MSSpectrum>(chromat, 1.0, 3.0, true);
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~SignalToNoiseOpenMS())
{
  delete ptr;
}
END_SECTION
START_SECTION(double getValueAtRT(double RT))
{
  static const double arr1[] = 
  {
    200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340,
    350, 360, 370, 380, 390, 400, 410, 420, 430, 440, 450, 460, 470, 480, 490,
    500, 510, 520, 530, 540, 550, 560, 570, 580, 590
  };
  std::vector<double> mz (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
  static const double arr2[] = 
  {
    5.4332, 5.6189, 4.3025, 4.5705, 5.4538, 9.7202, 8.805, 8.5391, 6.6257,
    5.809, 6.5518, 7.9273, 5.3875, 9.826, 5.139, 5.8588, 0.7806, 4.2054,
    9.9171, 4.0198, 1.1462, 5.1042, 7.8318, 4.8553, 6.691, 4.2377, 7.2344,
    4.0124, 3.8565, 6.2867, 1.0817, 8.2412, 5.0589, 7.0478, 5.9388, 1.2747,
    2.4228, 4.909, 6.856, 1.9665
  };
  std::vector<double> intensity (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

  MSSpectrum s;
  for (Size i = 0; i < mz.size(); i++)
  {
    Peak1D p;
    p.setMZ(mz[i]);
    p.setIntensity(intensity[i]);
    s.push_back(p);
  }
  SignalToNoiseOpenMS<MSSpectrum> ff(s, 200, 50, true);

  double value200 = 0.987854524;
  double value210 = 1.02162;
  double value220 = 0.782272686;
  double value590 = 0.35754546252164;

  // test values between the mz values
  TEST_REAL_SIMILAR ( ff.getValueAtRT(201), value200 )
  TEST_REAL_SIMILAR ( ff.getValueAtRT(211), value210 )
  TEST_REAL_SIMILAR ( ff.getValueAtRT(221), value220 )

  // test values exactly on the mz values
  TEST_REAL_SIMILAR ( ff.getValueAtRT(200), value200 )
  TEST_REAL_SIMILAR ( ff.getValueAtRT(210), value210 )

  // test values outside the range
  TEST_REAL_SIMILAR ( ff.getValueAtRT(100), value200 )
  TEST_REAL_SIMILAR ( ff.getValueAtRT(588), value590 )
  TEST_REAL_SIMILAR ( ff.getValueAtRT(590), value590 )
  TEST_REAL_SIMILAR ( ff.getValueAtRT(700), value590 )

}
END_SECTION


}
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

