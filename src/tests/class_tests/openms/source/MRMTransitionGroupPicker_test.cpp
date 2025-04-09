// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
///////////////////////////

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <boost/assign/std/vector.hpp>

using namespace OpenMS;
using namespace std;

typedef MSChromatogram RichPeakChromatogram;
// TODO also test the picker with the LightTransition interface
// typedef MRMTransitionGroup<RichPeakChromatogram, OpenSwath::LightTransition> MRMTransitionGroupType;
typedef MRMTransitionGroup<RichPeakChromatogram, ReactionMonitoringTransition> MRMTransitionGroupType;

void setup_transition_group(MRMTransitionGroupType & transition_group)
{
  // this is a simulated SRM experiment where the two traces are not sampled at
  // the exact same time points, thus a re-sampling is necessary before applying
  // the algorithm.
  // The MS1 trace intensity is a simple quadratic function.
  //

/*
 * Python code to create the MS1 trace :
 *
 
datapoints = [-100*(x-9)*(x-9)+9000 for x in range(18) ] 
sum(datapoints[3:10])
53900

*/

  static const double rtdata_1[] = {1474.34, 1477.11, 1479.88, 1482.64,
    1485.41, 1488.19, 1490.95, 1493.72, 1496.48, 1499.25, 1502.03, 1504.8,
    1507.56, 1510.33, 1513.09, 1515.87, 1518.64, 1521.42};
  static const double rtdata_2[] = {1473.55, 1476.31, 1479.08, 1481.84,
    1484.61, 1487.39, 1490.15, 1492.92, 1495.69, 1498.45, 1501.23, 1504,
    1506.76, 1509.53, 1512.29, 1515.07, 1517.84, 1520.62};

  static const double intdata_1[] = {3.26958, 3.74189, 3.31075, 86.1901,
    3.47528, 387.864, 13281  , 6375.84, 39852.6, 2.66726, 612.747, 3.34313,
    793.12 , 3.29156, 4.00586, 4.1591 , 3.23035, 3.90591};
  static const double intdata_2[] =  {3.44054 , 2142.31 , 3.58763 , 3076.97 ,
    6663.55 , 45681   , 157694  , 122844  , 86034.7 , 85391.1 , 15992.8 ,
    2293.94 , 6934.85 , 2735.18 , 459.413 , 3.93863 , 3.36564 , 3.44005};
  static const double ms1_intdata[] =  {900, 2600, 4100, 5400, 6500, 7400,
    8100, 8600, 8900, 9000, 8900, 8600, 8100, 7400, 6500, 5400, 4100, 2600};

  // Transition trace 1
  {
  ReactionMonitoringTransition transition;
  transition.setNativeID("1");
  RichPeakChromatogram chromatogram;
  for (int k = 0; k < 18; k++)
  {   
    ChromatogramPeak peak;
    peak.setRT(rtdata_1[k]);
    peak.setIntensity(intdata_1[k]);
    chromatogram.push_back(peak);
  } 
  chromatogram.setMetaValue("product_mz", 618.31);
  chromatogram.setNativeID("1");
  transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
  transition_group.addTransition(transition, transition.getNativeID());
  }

  // Transition trace 2
  {
  ReactionMonitoringTransition transition;
  transition.setNativeID("2");
  RichPeakChromatogram chromatogram;
  for (int k = 0; k < 18; k++)
  {   
    ChromatogramPeak peak;
    peak.setRT(rtdata_2[k]);
    peak.setIntensity(intdata_2[k]);
    chromatogram.push_back(peak);
  } 
  chromatogram.setMetaValue("product_mz", 619.31);
  chromatogram.setNativeID("2");
  transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
  transition_group.addTransition(transition, transition.getNativeID());
  }

  // MS1 trace
  {
  RichPeakChromatogram chromatogram;
  for (int k = 0; k < 18; k++)
  {   
    ChromatogramPeak peak;
    peak.setRT(rtdata_2[k] + 0.5); // shift the "MS1" retention time as well 
    peak.setIntensity(ms1_intdata[k]);
    chromatogram.push_back(peak);
  } 
  chromatogram.setNativeID("Precursor_i0");
  transition_group.addPrecursorChromatogram(chromatogram, "Precursor_i0");
  }
}

void setup_transition_group2(MRMTransitionGroupType & transition_group)
{
  // Simulated SRM experiment with a large peak in trace 1,
  // two smaller peaks in trace 2 that overlap with the peak in trace 1,
  // and no peaks in trace 3

  static const double rtdata_1[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
  static const double rtdata_2[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
  static const double rtdata_3[] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};

  static const double intdata_1[] = {3,3,3,3,3,3,3,3,4,6,8,12,14,15,16,15,14,12,8,6,4,3,3,3,3,3,3,3,3,3};
  static const double intdata_2[] = {3,3,4,6,8,10,11,11,10,8,6,4,4,6,8,10,11,11,10,8,6,4,3,3,3,3,3,3,3,3};
  static const double intdata_3[] =  {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};

  // Transition trace 1
  {
  ReactionMonitoringTransition transition;
  transition.setNativeID("1");
  RichPeakChromatogram chromatogram;
  for (int k = 0; k < 30; k++)
  {   
    ChromatogramPeak peak;
    peak.setRT(rtdata_1[k]);
    peak.setIntensity(intdata_1[k]);
    chromatogram.push_back(peak);
  } 
  chromatogram.setMetaValue("product_mz", 618.31);
  chromatogram.setNativeID("1");
  transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
  transition_group.addTransition(transition, transition.getNativeID());
  }

  // Transition trace 2
  {
  ReactionMonitoringTransition transition;
  transition.setNativeID("2");
  RichPeakChromatogram chromatogram;
  for (int k = 0; k < 30; k++)
  {   
    ChromatogramPeak peak;
    peak.setRT(rtdata_2[k]);
    peak.setIntensity(intdata_2[k]);
    chromatogram.push_back(peak);
  } 
  chromatogram.setMetaValue("product_mz", 619.31);
  chromatogram.setNativeID("2");
  transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
  transition_group.addTransition(transition, transition.getNativeID());
  }

  // Transition trace 3
  {
  ReactionMonitoringTransition transition;
  transition.setNativeID("3");
  RichPeakChromatogram chromatogram;
  for (int k = 0; k < 30; k++)
  {   
    ChromatogramPeak peak;
    peak.setRT(rtdata_3[k]);
    peak.setIntensity(intdata_3[k]);
    chromatogram.push_back(peak);
  } 
  chromatogram.setMetaValue("product_mz", 812.2);
  chromatogram.setNativeID("3");
  transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
  transition_group.addTransition(transition, transition.getNativeID());
  }
}

void setup_toy_chromatogram(RichPeakChromatogram & chromatogram)
{
  
  // Toy chromatogram
  // data is taken from raw LC-MS/MS data points acquired for L-Glutamate in RBCs
  std::vector<double> time={2.23095,2.239716667,2.248866667,2.25765,2.266416667,
    2.275566667,2.2847,2.293833333,2.304066667,2.315033333,2.325983333,2.336566667,
    2.3468,2.357016667,2.367283333,2.377183333,2.387083333,2.39735,2.40725,2.4175,
    2.4274,2.4373,2.44755,2.45745,2.4677,2.477966667,2.488216667,2.498516667,2.5084,
    2.5183,2.5282,2.538466667,2.548366667,2.558266667,2.568516667,2.578783333,
    2.588683333,2.59895,2.6092,2.619466667,2.630066667,2.64065,2.65125,2.662116667,
    2.672716667,2.6833,2.6939,2.7045,2.715083333,2.725683333,2.736266667,2.746866667,
    2.757833333,2.768416667,2.779016667,2.789616667,2.8002,2.810116667,2.820033333,
    2.830316667,2.840216667,2.849766667,2.859316667,2.868866667,2.878783333,2.888683333,
    2.898233333,2.907783333,2.916033333,2.924266667,2.93215,2.940383333,2.947933333,
    2.955816667,2.964066667,2.97195,2.979833333,2.987716667,2.995616667,3.003516667,
    3.011416667,3.01895,3.026833333,3.034366667,3.042266667,3.0498,3.05735,3.065233333,
    3.073133333,3.080666667,3.0882,3.095733333,3.103633333,3.111533333,3.119066667,
    3.126966667,3.134866667,3.14275,3.15065,3.15855,3.166433333,3.174333333,3.182233333,
    3.190133333,3.198016667,3.205916667,3.213166667};
  std::vector<double> intensity={1447,2139,1699,755,1258,1070,944,1258,1573,1636,
    1762,1447,1133,1321,1762,1133,1447,2391,692,1636,2957,1321,1573,1196,1258,881,
    1384,2076,1133,1699,1384,692,1636,1133,1573,1825,1510,2391,4342,10382,17618,
    51093,153970,368094,632114,869730,962547,966489,845055,558746,417676,270942,
    184865,101619,59776,44863,31587,24036,20450,20324,11074,9879,10508,7928,7110,
    6733,6481,5726,6921,6670,5537,4971,4719,4782,5097,5789,4279,5411,4530,3524,
    2139,3335,3083,4342,4279,3083,3649,4216,4216,3964,2957,2202,2391,2643,3524,
    2328,2202,3649,2706,3020,3335,2580,2328,2894,3146,2769,2517};

  for (size_t k = 0; k < time.size(); k++)
  {   
    ChromatogramPeak peak;
    peak.setRT(time[k]);
    peak.setIntensity(intensity[k]);
    chromatogram.push_back(peak);
  } 
}

START_TEST(MRMTransitionGroupPicker, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMTransitionGroupPicker* ptr = nullptr;
MRMTransitionGroupPicker* nullPointer = nullptr;

START_SECTION(MRMTransitionGroupPicker())
{
  ptr = new MRMTransitionGroupPicker();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMTransitionGroupPicker())
{
  delete ptr;
}
END_SECTION

START_SECTION((template < typename SpectrumT, typename TransitionT > void pickTransitionGroup(MRMTransitionGroup< SpectrumT, TransitionT > &transition_group)))
{ 
  { // transition group 1
    MRMTransitionGroupType transition_group;
    setup_transition_group(transition_group);

    MRMTransitionGroupPicker trgroup_picker;
    Param picker_param = trgroup_picker.getDefaults();
    picker_param.setValue("PeakPickerChromatogram:method", "legacy"); // old parameters
    picker_param.setValue("PeakPickerChromatogram:peak_width", 40.0); // old parameters
    trgroup_picker.setParameters(picker_param);
    trgroup_picker.pickTransitionGroup(transition_group);

    TEST_EQUAL(transition_group.getFeatures().size(), 1)
    MRMFeature mrmfeature = transition_group.getFeatures()[0];
    TEST_REAL_SIMILAR(mrmfeature.getRT(), 1492.83060)
    TEST_REAL_SIMILAR(mrmfeature.getIntensity(), 567375)
    TEST_REAL_SIMILAR(mrmfeature.getMetaValue("leftWidth"), 1481.84)
    TEST_REAL_SIMILAR(mrmfeature.getMetaValue("rightWidth"), 1501.23)

    // test the number of hull points (should be equal)
    TEST_EQUAL( mrmfeature.getFeature("1").getConvexHulls()[0].getHullPoints().size(), 7);
    TEST_EQUAL( mrmfeature.getFeature("2").getConvexHulls()[0].getHullPoints().size(), 7);

    // the intensity of the hull points should not have changed

    // Check Feature 2
    ConvexHull2D::PointArrayType data1_points = mrmfeature.getFeature("2").getConvexHulls()[0].getHullPoints();
    double sum = 0.0;
    for (ConvexHull2D::PointArrayType::iterator it = data1_points.begin(); it != data1_points.end(); it++)
    {
      sum += it->getY();
    }
    TEST_REAL_SIMILAR(sum, 507385.32);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity(), 507385.32);

    // Check Feature 1
    ConvexHull2D::PointArrayType data2_points = mrmfeature.getFeature("1").getConvexHulls()[0].getHullPoints();
    sum = 0.0;
    for (ConvexHull2D::PointArrayType::iterator it = data2_points.begin(); it != data2_points.end(); it++)
    {
      sum += it->getY();
    }
    TEST_REAL_SIMILAR(sum, 59989.8287208466);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 59989.8287208466);

    // Also check the MS1
    std::vector<String> result;
    mrmfeature.getPrecursorFeatureIDs(result);
    TEST_EQUAL(result.size(), 1);
    TEST_EQUAL(result[0], "Precursor_i0");
    data2_points = mrmfeature.getPrecursorFeature("Precursor_i0").getConvexHulls()[0].getHullPoints();
    sum = 0.0;
    for (ConvexHull2D::PointArrayType::iterator it = data2_points.begin(); it != data2_points.end(); it++)
    {
      sum += it->getY();
    }
    // Part of the signal gets lost due to the re-sampling since the MS1 sampling
    // positions are not at the same place as the MS2 sampling positions
    double resampling_loss = 875.9514;
    TEST_REAL_SIMILAR(sum, 53900 - resampling_loss);
    TEST_REAL_SIMILAR(mrmfeature.getPrecursorFeature("Precursor_i0").getIntensity(),
    // mrmfeature.getMS1Feature().getIntensity(), 53900 - resampling_loss);
    53900 - resampling_loss);
  }

  { // transition group 2 (Consensus)
    MRMTransitionGroupType transition_group;
    setup_transition_group2(transition_group);

    MRMTransitionGroupPicker trgroup_picker;
    Param picker_param = trgroup_picker.getDefaults();
    picker_param.setValue("resample_boundary", 1.0);
    picker_param.setValue("PeakPickerChromatogram:gauss_width", 10.0);
    picker_param.setValue("PeakPickerChromatogram:peak_width", -1.0);
    picker_param.setValue("PeakPickerChromatogram:signal_to_noise", 1.0);
    picker_param.setValue("PeakPickerChromatogram:method", "corrected");
    picker_param.setValue("use_consensus", "true");
    trgroup_picker.setParameters(picker_param);
    trgroup_picker.pickTransitionGroup(transition_group);

    TEST_EQUAL(transition_group.getFeatures().size(), 1);
    MRMFeature mrmfeature = transition_group.getFeatures()[0];
    TEST_REAL_SIMILAR(mrmfeature.getRT(), 14);
    TEST_REAL_SIMILAR(mrmfeature.getMetaValue("leftWidth"), 7);
    TEST_REAL_SIMILAR(mrmfeature.getMetaValue("rightWidth"), 21);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 140);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity(), 117);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("3").getIntensity(), 45);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getMetaValue("peak_apex_int"), 16);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getMetaValue("peak_apex_int"), 11);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("3").getMetaValue("peak_apex_int"), 3);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getMetaValue("peak_apex_position"), 14);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getMetaValue("peak_apex_position"), 7);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("3").getMetaValue("peak_apex_position"), 7);
  }

  { // transition group 2 (no consensus)
    MRMTransitionGroupType transition_group;
    setup_transition_group2(transition_group);

    MRMTransitionGroupPicker trgroup_picker;
    Param picker_param = trgroup_picker.getDefaults();
    picker_param.setValue("resample_boundary", 1.0);
    picker_param.setValue("PeakPickerChromatogram:gauss_width", 10.0);
    picker_param.setValue("PeakPickerChromatogram:peak_width", -1.0);
    picker_param.setValue("PeakPickerChromatogram:signal_to_noise", 1.0);
    picker_param.setValue("PeakPickerChromatogram:method", "corrected");
    picker_param.setValue("use_consensus", "false");
    trgroup_picker.setParameters(picker_param);
    trgroup_picker.pickTransitionGroup(transition_group);

    TEST_EQUAL(transition_group.getFeatures().size(), 2);
    MRMFeature mrmfeature = transition_group.getFeatures()[0];
    TEST_REAL_SIMILAR(mrmfeature.getRT(), 14);
    TEST_REAL_SIMILAR(mrmfeature.getMetaValue("leftWidth"), 7);
    TEST_REAL_SIMILAR(mrmfeature.getMetaValue("rightWidth"), 21);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 140);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity(), 58); // consensus = 117
    TEST_REAL_SIMILAR(mrmfeature.getFeature("3").getIntensity(), 45);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getMetaValue("peak_apex_int"), 16);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getMetaValue("peak_apex_int"), 11);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("3").getMetaValue("peak_apex_int"), 3);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getMetaValue("peak_apex_position"), 14);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getMetaValue("peak_apex_position"), 16); // consensus = 7
    TEST_REAL_SIMILAR(mrmfeature.getFeature("3").getMetaValue("peak_apex_position"), 7);
  }

  { // transition group 1 -- only quantifying
    MRMTransitionGroupPicker trgroup_picker;
    Param picker_param = trgroup_picker.getDefaults();
    picker_param.setValue("PeakPickerChromatogram:method", "legacy"); // old parameters
    picker_param.setValue("PeakPickerChromatogram:peak_width", 40.0); // old parameters
    trgroup_picker.setParameters(picker_param);

    {
      MRMTransitionGroupType transition_group;
      setup_transition_group(transition_group);
      trgroup_picker.pickTransitionGroup(transition_group);

      TEST_EQUAL(transition_group.getFeatures().size(), 1)
      MRMFeature mrmfeature = transition_group.getFeatures()[0];
      TEST_REAL_SIMILAR(mrmfeature.getIntensity(), 567375)
      TEST_REAL_SIMILAR(mrmfeature.getRT(), 1492.83060);
    }

    {
      MRMTransitionGroupType transition_group;
      setup_transition_group(transition_group);
      transition_group.getTransitionsMuteable()[0].setQuantifyingTransition(false);
      trgroup_picker.pickTransitionGroup(transition_group);

      TEST_EQUAL(transition_group.getFeatures().size(), 1)
      MRMFeature mrmfeature = transition_group.getFeatures()[0];
      TEST_REAL_SIMILAR(mrmfeature.getIntensity(), 507385)
      TEST_REAL_SIMILAR(mrmfeature.getRT(), 1492.83060);
    }

    {
      MRMTransitionGroupType transition_group;
      setup_transition_group(transition_group);
      transition_group.getTransitionsMuteable()[1].setQuantifyingTransition(false);
      trgroup_picker.pickTransitionGroup(transition_group);

      TEST_EQUAL(transition_group.getFeatures().size(), 1)
      MRMFeature mrmfeature = transition_group.getFeatures()[0];
      TEST_REAL_SIMILAR(mrmfeature.getIntensity(), 59989.8)
      TEST_REAL_SIMILAR(mrmfeature.getRT(), 1492.83060);
    }

    {
      MRMTransitionGroupType transition_group;
      setup_transition_group(transition_group);
      transition_group.getTransitionsMuteable()[0].setQuantifyingTransition(false);
      transition_group.getTransitionsMuteable()[1].setQuantifyingTransition(false);
      trgroup_picker.pickTransitionGroup(transition_group);

      TEST_EQUAL(transition_group.getFeatures().size(), 0) // no results, since zero intensity transitions get removed
    }
  }

  { // transition group 1 -- only detecting
    MRMTransitionGroupPicker trgroup_picker;
    Param picker_param = trgroup_picker.getDefaults();
    picker_param.setValue("resample_boundary", 1.0);
    picker_param.setValue("PeakPickerChromatogram:gauss_width", 10.0);
    picker_param.setValue("PeakPickerChromatogram:peak_width", -1.0);
    picker_param.setValue("PeakPickerChromatogram:signal_to_noise", 1.0);
    picker_param.setValue("PeakPickerChromatogram:method", "corrected");
    picker_param.setValue("use_consensus", "false");
    trgroup_picker.setParameters(picker_param);

    {
      MRMTransitionGroupType transition_group;
      setup_transition_group2(transition_group);
      trgroup_picker.pickTransitionGroup(transition_group);

      TEST_EQUAL(transition_group.getFeatures().size(), 2);
      MRMFeature mrmfeature = transition_group.getFeatures()[0];
      TEST_REAL_SIMILAR(mrmfeature.getIntensity(), 243)
      TEST_REAL_SIMILAR(mrmfeature.getRT(), 14.0);

      TEST_REAL_SIMILAR(mrmfeature.getRT(), 14);
      TEST_REAL_SIMILAR(mrmfeature.getMetaValue("leftWidth"), 7);
      TEST_REAL_SIMILAR(mrmfeature.getMetaValue("rightWidth"), 21);
      TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 140);
      TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity(), 58); // consensus = 117
      TEST_REAL_SIMILAR(mrmfeature.getFeature("3").getIntensity(), 45);
      TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getMetaValue("peak_apex_int"), 16);
      TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getMetaValue("peak_apex_int"), 11);
      TEST_REAL_SIMILAR(mrmfeature.getFeature("3").getMetaValue("peak_apex_int"), 3);
      TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getMetaValue("peak_apex_position"), 14);
      TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getMetaValue("peak_apex_position"), 16); // consensus = 7
      TEST_REAL_SIMILAR(mrmfeature.getFeature("3").getMetaValue("peak_apex_position"), 7);
    }

    {
      // cannot have non-consensus data where we wont use all transitions
      MRMTransitionGroupType transition_group;
      setup_transition_group2(transition_group);
      transition_group.getTransitionsMuteable()[0].setDetectingTransition(false);
      TEST_EXCEPTION(Exception::IllegalArgument, trgroup_picker.pickTransitionGroup(transition_group))
      TEST_EQUAL(transition_group.getFeatures().size(), 0)
    }

  }
}
END_SECTION

START_SECTION((template <typename SpectrumT, typename TransitionT> MRMFeature createMRMFeature(MRMTransitionGroup<SpectrumT, TransitionT>& transition_group, std::vector<SpectrumT>& picked_chroms, std::vector<SpectrumT>& smoothed_chroms, int& chr_idx, int& peak_idx)))
{
  MRMTransitionGroupType transition_group;
  setup_transition_group(transition_group);
  std::vector< RichPeakChromatogram > picked_chroms;
  std::vector< RichPeakChromatogram > smoothed_chroms;

  double left_start = 1481.840;
  double right_end = 1512.290;

  // do "peakpicking", create one peak
  for(Size k = 0; k < transition_group.getChromatograms().size(); k++)
  {
    RichPeakChromatogram picked_chrom;
    ChromatogramPeak peak;
    peak.setRT(1490);
    peak.setIntensity(170);
    picked_chrom.push_back(peak);

    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(PeakPickerChromatogram::SIZE_OF_FLOATINDICES);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].setName("rightWidth");
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE].push_back(1000.0);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(left_start);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(right_end);
    picked_chrom.setNativeID(transition_group.getChromatograms()[k].getNativeID());

    picked_chroms.push_back(picked_chrom);
  }

  // create the corresponding first mrm feature
  int chr_idx = 1, peak_idx = 0;

  {
    MRMTransitionGroupPicker picker;
    Param picker_param = picker.getDefaults();
    picker_param.setValue("PeakPickerChromatogram::method", "legacy"); // old parameters
    picker_param.setValue("PeakPickerChromatogram::peak_width", 40.0); // old parameters
    picker.setParameters(picker_param);

    MRMFeature mrmfeature = picker.createMRMFeature(transition_group, picked_chroms, smoothed_chroms, chr_idx, peak_idx);
    TEST_REAL_SIMILAR(mrmfeature.getRT(), 1490.0)

    // test the number of hull points (should be equal)
    TEST_EQUAL( mrmfeature.getFeature("1").getConvexHulls()[0].getHullPoints().size(), 12);
    TEST_EQUAL( mrmfeature.getFeature("2").getConvexHulls()[0].getHullPoints().size(), 12);

    // the intensity of the hull points should not have changed
    ConvexHull2D::PointArrayType data1_points = mrmfeature.getFeature("2").getConvexHulls()[0].getHullPoints();
    double sum = 0.0;
    for (ConvexHull2D::PointArrayType::iterator it = data1_points.begin(); it != data1_points.end(); it++)
    {
      sum += it->getY();
    }
    TEST_REAL_SIMILAR(sum, 535801.503);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity(), 535801.503);
    TEST_REAL_SIMILAR((double)mrmfeature.getFeature("2").getMetaValue("peak_apex_int"), 157694);
    ConvexHull2D::PointArrayType data2_points = mrmfeature.getFeature("1").getConvexHulls()[0].getHullPoints();
    sum = 0.0;
    for (ConvexHull2D::PointArrayType::iterator it = data2_points.begin(); it != data2_points.end(); it++)
    {
      sum += it->getY();
    }
    TEST_REAL_SIMILAR(sum, 61405.95106);
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 61405.95106);
    TEST_REAL_SIMILAR((double)mrmfeature.getFeature("1").getMetaValue("peak_apex_int"), 30286.9130513267);

    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity()/mrmfeature.getFeature("1").getIntensity(), 8.72556) // ratio should be stable

    // feature dimension
    TEST_EQUAL(mrmfeature.getRT(), 1490.0)
    TEST_REAL_SIMILAR( mrmfeature.getMetaValue("leftWidth"), left_start);
    TEST_REAL_SIMILAR( mrmfeature.getMetaValue("rightWidth"), right_end);
  }

  {
    MRMTransitionGroupPicker picker;
    Param picker_param = picker.getDefaults();
    picker_param.setValue("PeakPickerChromatogram::method", "legacy"); // old parameters
    picker_param.setValue("PeakPickerChromatogram::peak_width", 40.0); // old parameters
    picker_param.setValue("background_subtraction", "original");
    picker.setParameters(picker_param);

    MRMFeature mrmfeature = picker.createMRMFeature(transition_group, picked_chroms, smoothed_chroms, chr_idx, peak_idx);

    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity(), 514583); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR((double)mrmfeature.getFeature("2").getMetaValue("peak_apex_int"), 155925.8085); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 60913.8); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR((double)mrmfeature.getFeature("2").getMetaValue("peak_apex_int"), 155925.8085); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity()/mrmfeature.getFeature("1").getIntensity(), 8.44773) // ratio should be stable
  }

  // integration with Simpon's rule (with background subtraction)
  {
    MRMTransitionGroupPicker picker;
    Param picker_param = picker.getDefaults();
    picker_param.setValue("PeakPickerChromatogram::method", "legacy"); // old parameters
    picker_param.setValue("PeakPickerChromatogram::peak_width", 40.0); // old parameters
    picker_param.setValue("background_subtraction", "exact");
    picker_param.setValue("PeakIntegrator:integration_type", "simpson");
    picker.setParameters(picker_param);

    MRMFeature mrmfeature = picker.createMRMFeature(transition_group, picked_chroms, smoothed_chroms, chr_idx, peak_idx);

    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity(), 1.42194e+06); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR((double)mrmfeature.getFeature("2").getMetaValue("peak_apex_int"), 155331.37806798); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 168685); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR((double)mrmfeature.getFeature("1").getMetaValue("peak_apex_int"), 30251.2414481247); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity()/mrmfeature.getFeature("1").getIntensity(), 8.42957) // ratio should be stable
  }

  // background subtraction with peak integrator
  {
    MRMTransitionGroupPicker picker;
    Param picker_param = picker.getDefaults();
    picker_param.setValue("PeakPickerChromatogram::method", "legacy"); // old parameters
    picker_param.setValue("PeakPickerChromatogram::peak_width", 40.0); // old parameters
    picker_param.setValue("background_subtraction", "exact");
    picker.setParameters(picker_param);

    MRMFeature mrmfeature = picker.createMRMFeature(transition_group, picked_chroms, smoothed_chroms, chr_idx, peak_idx);

    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity(), 514590); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR((double)mrmfeature.getFeature("2").getMetaValue("peak_apex_int"), 155331.37806798); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 61009.7); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR((double)mrmfeature.getFeature("1").getMetaValue("peak_apex_int"), 30251.2414481247); // slightly reduced value due to background subtraction
    TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity()/mrmfeature.getFeature("1").getIntensity(), 8.43456040596823) // ratio should be stable
  }
}
END_SECTION

///////////////////////////////////////////////////////////////////////////
/// Private methods
///////////////////////////////////////////////////////////////////////////

START_SECTION(( void findLargestPeak(std::vector<RichPeakChromatogram> & picked_chroms, int & chr_idx, int & peak_idx) ))
{
  MRMTransitionGroupType transition_group;
  std::vector< RichPeakChromatogram > picked_chroms;

  // do peakpicking, e.g. find a peak at 3120 RT / 170 intensity in all the spectra
  for(Size k = 0; k < 3; k++)
  {
    RichPeakChromatogram picked_chrom;
    ChromatogramPeak peak;
    peak.setRT(3120);
    peak.setIntensity(100+k);
    picked_chrom.push_back(peak);
    peak.setRT(4120);
    peak.setIntensity(200+k);
    picked_chrom.push_back(peak);

    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(PeakPickerChromatogram::SIZE_OF_FLOATINDICES);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].setName("rightWidth");
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_ABUNDANCE].push_back(1000.0);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(3100.0);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(3140.0);

    picked_chroms.push_back(picked_chrom);
  }

  MRMTransitionGroupPicker picker;
  int chr_idx = -1, peak_idx = -1;
  picker.findLargestPeak(picked_chroms, chr_idx, peak_idx);

  TEST_EQUAL(chr_idx, 2);
  TEST_EQUAL(peak_idx, 1);
}
END_SECTION

START_SECTION(void findWidestPeakIndices(const std::vector<MSChromatogram>& picked_chroms, Int& chrom_idx, Int& point_idx) const)
{
  std::vector<MSChromatogram> chromatograms;
  MSChromatogram c;
  c.push_back(ChromatogramPeak(110.0, 2000.0));
  c.push_back(ChromatogramPeak(150.0, 6000.0));
  c.push_back(ChromatogramPeak(190.0, 2000.0));
  c.getFloatDataArrays().resize(PeakPickerChromatogram::SIZE_OF_FLOATINDICES);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(100.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(120.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(120.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(180.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(180.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(200.0);
  chromatograms.push_back(c); // chromatogram containing a peak of highest intensity (should be skipped in favor of widest peak)

  c.clear(true);
  c.push_back(ChromatogramPeak(150.0, 5500.0)); // lower global intensity, if compared to the previous chromatogram
  c.push_back(ChromatogramPeak(190.0, 2000.0));
  c.getFloatDataArrays().resize(PeakPickerChromatogram::SIZE_OF_FLOATINDICES);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(100.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(180.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(180.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(200.0);
  chromatograms.push_back(c); // chromatogram containing the widest peak (this should be chosen)

  c.clear(true);
  c.push_back(ChromatogramPeak(110.0, 2000.0));
  c.push_back(ChromatogramPeak(150.0, 7000.0));
  c.push_back(ChromatogramPeak(190.0, 2000.0));
  c.getFloatDataArrays().resize(PeakPickerChromatogram::SIZE_OF_FLOATINDICES);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(105.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(115.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(125.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(175.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(185.0);
  c.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(195.0);
  chromatograms.push_back(c); // just another chromatogram (it won't be chosen, it contains short peaks)

  MRMTransitionGroupPicker picker;
  Int chr_idx{-1}, peak_idx{-1};
  picker.findWidestPeakIndices(chromatograms, chr_idx, peak_idx);
  TEST_EQUAL(chr_idx, 1);
  TEST_EQUAL(peak_idx, 0); // the point [0] (first) is the apex of the widest peak within the chosen chromatogram
  TEST_REAL_SIMILAR(chromatograms[chr_idx][peak_idx].getRT(), 150.0)
  TEST_REAL_SIMILAR(chromatograms[chr_idx][peak_idx].getIntensity(), 5500.0)
  TEST_REAL_SIMILAR(chromatograms[chr_idx].getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][peak_idx], 100.0)
  TEST_REAL_SIMILAR(chromatograms[chr_idx].getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][peak_idx], 180.0)
}
END_SECTION

START_SECTION((template < typename SpectrumT > void remove_overlapping_features(std::vector< SpectrumT > &picked_chroms, double best_left, double best_right)))
{
  MRMTransitionGroupType transition_group;
  std::vector< RichPeakChromatogram > picked_chroms;
  MRMTransitionGroupPicker picker;
  double default_intensity = 170;

  // create 3 peaks, at 3120, 3090 and 3060 which are all overlapping
  {
    RichPeakChromatogram picked_chrom;
    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(PeakPickerChromatogram::SIZE_OF_FLOATINDICES);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].setName("rightWidth");

    {
    ChromatogramPeak peak;
    peak.setRT(3120);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(3100.0);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(3140.0);
    }

    {
    ChromatogramPeak peak;
    peak.setRT(3090);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(3070.0);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(3120.0);
    }

    {
    ChromatogramPeak peak;
    peak.setRT(3060);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(3050.0);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(3090.0);
    }
    
    picked_chroms.push_back(picked_chrom);
  }

  // create 2 peaks, at 3120 and 3060 which are not overlapping
  {
    RichPeakChromatogram picked_chrom;
    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(PeakPickerChromatogram::SIZE_OF_FLOATINDICES);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].setName("rightWidth");

    {
    ChromatogramPeak peak;
    peak.setRT(3120);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(3100.0);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(3140.0);
    }

    {
    ChromatogramPeak peak;
    peak.setRT(3060);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER].push_back(3050.0);
    picked_chrom.getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER].push_back(3090.0);
    }
    
    picked_chroms.push_back(picked_chrom);
  }

  std::vector< RichPeakChromatogram > picked_chroms_orig = picked_chroms;

  // First we look at the rightmost peak which should include the first two
  // peaks in the first chromatogram and the first peak in the second
  //chromatogram
  double best_left = 3100; double best_right = 3140; 
  picked_chroms = picked_chroms_orig;
  picker.remove_overlapping_features(picked_chroms, best_left, best_right);

  TEST_REAL_SIMILAR(picked_chroms[0][0].getIntensity(), 0.0)
  TEST_REAL_SIMILAR(picked_chroms[0][1].getIntensity(), 0.0)
  TEST_REAL_SIMILAR(picked_chroms[0][2].getIntensity(), default_intensity)
  TEST_REAL_SIMILAR(picked_chroms[1][0].getIntensity(), 0.0)
  TEST_REAL_SIMILAR(picked_chroms[1][1].getIntensity(), default_intensity)

  // Second we look at the middle peak which should include all peaks
  best_left = 3070; best_right = 3120; 
  picked_chroms = picked_chroms_orig;
  picker.remove_overlapping_features(picked_chroms, best_left, best_right);

  TEST_REAL_SIMILAR(picked_chroms[0][0].getIntensity(), 0.0)
  TEST_REAL_SIMILAR(picked_chroms[0][1].getIntensity(), 0.0)
  TEST_REAL_SIMILAR(picked_chroms[0][2].getIntensity(), 0.0);
  TEST_REAL_SIMILAR(picked_chroms[1][0].getIntensity(), 0.0)
  TEST_REAL_SIMILAR(picked_chroms[1][1].getIntensity(), 0.0);

  // Last we look at the leftmost peak which should include all peaks
  best_left = 3050; best_right = 3090; 
  picked_chroms = picked_chroms_orig;
  picker.remove_overlapping_features(picked_chroms, best_left, best_right);

  TEST_REAL_SIMILAR(picked_chroms[0][0].getIntensity(), default_intensity)
  TEST_REAL_SIMILAR(picked_chroms[0][1].getIntensity(), 0.0)
  TEST_REAL_SIMILAR(picked_chroms[0][2].getIntensity(), 0.0);
  TEST_REAL_SIMILAR(picked_chroms[1][0].getIntensity(), default_intensity)
  TEST_REAL_SIMILAR(picked_chroms[1][1].getIntensity(), 0.0);

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
