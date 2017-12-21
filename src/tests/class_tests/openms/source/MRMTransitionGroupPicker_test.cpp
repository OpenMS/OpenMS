// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    1485.41, 1488.19, 1490.95, 1493.72, 1496.48, 1499.25, 1502.03, 1504.8 ,
    1507.56, 1510.33, 1513.09, 1515.87, 1518.64, 1521.42};
  static const double rtdata_2[] = {1473.55, 1476.31, 1479.08, 1481.84,
    1484.61, 1487.39, 1490.15, 1492.92, 1495.69, 1498.45, 1501.23, 1504   ,
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
    peak.setMZ(rtdata_1[k]);
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
    peak.setMZ(rtdata_2[k]);
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
    peak.setMZ(rtdata_2[k] + 0.5); // shift the "MS1" retention time as well 
    peak.setIntensity(ms1_intdata[k]);
    chromatogram.push_back(peak);
  } 
  chromatogram.setNativeID("Precursor_i0");
  transition_group.addPrecursorChromatogram(chromatogram, "Precursor_i0");
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
    peak.setMZ(time[k]);
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
  MRMTransitionGroupType transition_group;
  setup_transition_group(transition_group);

  MRMTransitionGroupPicker trgroup_picker;
  Param picker_param = trgroup_picker.getDefaults();
  picker_param.setValue("PeakPickerMRM:method", "legacy"); // old parameters
  picker_param.setValue("PeakPickerMRM:peak_width", 40.0); // old parameters
  trgroup_picker.setParameters(picker_param);
  trgroup_picker.pickTransitionGroup(transition_group);

  TEST_EQUAL(transition_group.getFeatures().size(), 1)
  MRMFeature mrmfeature = transition_group.getFeatures()[0];
  TEST_REAL_SIMILAR(mrmfeature.getRT(), 1492.83060);
  TEST_REAL_SIMILAR( mrmfeature.getMetaValue("leftWidth"), 1481.84);
  TEST_REAL_SIMILAR( mrmfeature.getMetaValue("rightWidth"), 1501.23);

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
    peak.setMZ(1490);
    peak.setIntensity(170);
    picked_chrom.push_back(peak);

    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(3);
    picked_chrom.getFloatDataArrays()[0].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[2].setName("rightWidth");
    picked_chrom.getFloatDataArrays()[0].push_back(1000.0);
    picked_chrom.getFloatDataArrays()[1].push_back(left_start);
    picked_chrom.getFloatDataArrays()[2].push_back(right_end);
    picked_chrom.setNativeID(transition_group.getChromatograms()[k].getNativeID());

    picked_chroms.push_back(picked_chrom);
  }

  // create the corresponding first mrm feature
  int chr_idx = 1, peak_idx = 0;
  MRMTransitionGroupPicker picker;

  Param picker_param = picker.getDefaults();
  picker_param.setValue("PeakPickerMRM::method", "legacy"); // old parameters
  picker_param.setValue("PeakPickerMRM::peak_width", 40.0); // old parameters
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
  ConvexHull2D::PointArrayType data2_points = mrmfeature.getFeature("1").getConvexHulls()[0].getHullPoints();
  sum = 0.0;
  for (ConvexHull2D::PointArrayType::iterator it = data2_points.begin(); it != data2_points.end(); it++)
  {
    sum += it->getY();
  }
  TEST_REAL_SIMILAR(sum, 61405.95106);
  TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 61405.95106);

  // feature dimension
  TEST_EQUAL(mrmfeature.getRT(), 1490.0)
  TEST_REAL_SIMILAR( mrmfeature.getMetaValue("leftWidth"), left_start);
  TEST_REAL_SIMILAR( mrmfeature.getMetaValue("rightWidth"), right_end);
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
    peak.setMZ(3120);
    peak.setIntensity(100+k);
    picked_chrom.push_back(peak);
    peak.setMZ(4120);
    peak.setIntensity(200+k);
    picked_chrom.push_back(peak);

    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(3);
    picked_chrom.getFloatDataArrays()[0].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[2].setName("rightWidth");
    picked_chrom.getFloatDataArrays()[0].push_back(1000.0);
    picked_chrom.getFloatDataArrays()[1].push_back(3100.0);
    picked_chrom.getFloatDataArrays()[2].push_back(3140.0);

    picked_chroms.push_back(picked_chrom);
  }

  MRMTransitionGroupPicker picker;
  int chr_idx = -1, peak_idx = -1;
  picker.findLargestPeak(picked_chroms, chr_idx, peak_idx);

  TEST_EQUAL(chr_idx, 2);
  TEST_EQUAL(peak_idx, 1);
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
    picked_chrom.getFloatDataArrays().resize(3);
    picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[2].setName("rightWidth");

    {
    ChromatogramPeak peak;
    peak.setMZ(3120);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[1].push_back(3100.0);
    picked_chrom.getFloatDataArrays()[2].push_back(3140.0);
    }

    {
    ChromatogramPeak peak;
    peak.setMZ(3090);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[1].push_back(3070.0);
    picked_chrom.getFloatDataArrays()[2].push_back(3120.0);
    }

    {
    ChromatogramPeak peak;
    peak.setMZ(3060);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[1].push_back(3050.0);
    picked_chrom.getFloatDataArrays()[2].push_back(3090.0);
    }
    
    picked_chroms.push_back(picked_chrom);
  }

  // create 2 peaks, at 3120 and 3060 which are not overlapping
  {
    RichPeakChromatogram picked_chrom;
    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(3);
    picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[2].setName("rightWidth");

    {
    ChromatogramPeak peak;
    peak.setMZ(3120);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[1].push_back(3100.0);
    picked_chrom.getFloatDataArrays()[2].push_back(3140.0);
    }

    {
    ChromatogramPeak peak;
    peak.setMZ(3060);
    peak.setIntensity(default_intensity);
    picked_chrom.push_back(peak);
    picked_chrom.getFloatDataArrays()[1].push_back(3050.0);
    picked_chrom.getFloatDataArrays()[2].push_back(3090.0);
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

START_SECTION(( void calculateBgEstimationAverage_(const MSChromatogram& chromatogram,
  double best_left, double best_right, double & background, double & avg_noise_level) ))
{

  RichPeakChromatogram chromatogram;
  setup_toy_chromatogram(chromatogram);

  // Features
  double best_left = 2.477966667;
  double best_right = 3.01895;
  double background, noise_level;

  // Correct the background
  MRMTransitionGroupPicker picker;

  picker.calculateBgEstimationAverage_(chromatogram, 
    best_left, best_right, background,
    noise_level);

  TEST_REAL_SIMILAR(background, 125076);
  TEST_REAL_SIMILAR(noise_level, 2233.5);
}
END_SECTION

START_SECTION(( void calculateBgEstimationExact_(const MSChromatogram& chromatogram,
  double best_left, double best_right, double peak_height, double & background, double & avg_noise_level) ))
{

  RichPeakChromatogram chromatogram;
  setup_toy_chromatogram(chromatogram);

  // Features
  double best_left = 2.477966667;
  double best_right = 3.01895;
  double peak_height = 966489;
  double background, noise_level;

  // Correct the background
  MRMTransitionGroupPicker picker;

  picker.calculateBgEstimationExact_(chromatogram, 
    best_left, best_right, 
    peak_height, background,
    noise_level);

  TEST_REAL_SIMILAR(background, 123446.661339019);
  TEST_REAL_SIMILAR(noise_level, 1908.596906);
}
END_SECTION

START_SECTION(( void calculatePeakApexInt_(const MSChromatogram& chromatogram,
  double best_left, double best_right, 
  ConvexHull2D::PointArrayType & hull_points,
  double & intensity_sum, 
  double & rt_sum,
  double & peak_apex_int,
  double & peak_apex_rt) ))
{

  RichPeakChromatogram chromatogram;
  setup_toy_chromatogram(chromatogram);

  // Features
  double best_left = 2.477966667;
  double best_right = 3.01895;
  double peak_apex = 2.7045;

  // Calculate peak apex
  MRMTransitionGroupPicker picker;  

  ConvexHull2D::PointArrayType hull_points;
  double intensity_integral(0.0), intensity_sum(0.0), rt_sum(0.0);
  double peak_apex_int = -1;

  picker.calculatePeakApexInt_(chromatogram,
    best_left,best_right,hull_points,
    intensity_sum,
    intensity_integral,
    rt_sum,
    peak_apex_int,
    peak_apex);

  TEST_REAL_SIMILAR(intensity_sum, 6764562);
  TEST_REAL_SIMILAR(intensity_integral, 71540.2082038256);
  TEST_REAL_SIMILAR(rt_sum, 151.890633338);
  TEST_REAL_SIMILAR(peak_apex_int, 966489);
}
END_SECTION

START_SECTION(( void calculatePeakShapeMetrics_(const MSChromatogram& chromatogram, 
  double best_left, double best_right, 
  double peak_height, double peak_apex, double avg_noise_level,
  PeakShapeMetrics_ & peakShapeMetrics) ))
{
  
  RichPeakChromatogram chromatogram;
  setup_toy_chromatogram(chromatogram);

  // Features
  double best_left = 2.477966667;
  double best_right = 3.01895;
  double peak_height = 965356;
  double peak_apex = 2.7045;
  double avg_noise_level = 723.5;

  // Calculate the QCs
  MRMTransitionGroupPicker picker;
  MRMTransitionGroupPicker::PeakShapeMetrics_ peakShapeMetrics;

  picker.calculatePeakShapeMetrics_(chromatogram, 
    best_left, best_right, 
    peak_height, peak_apex, avg_noise_level,
    peakShapeMetrics);

  TEST_REAL_SIMILAR(peakShapeMetrics.width_at_5,0.27924231787346);
  TEST_REAL_SIMILAR(peakShapeMetrics.width_at_10,0.135162753574054);
  TEST_REAL_SIMILAR(peakShapeMetrics.width_at_50,0.0596533918928945);
  TEST_REAL_SIMILAR(peakShapeMetrics.start_time_at_10,2.63202095937465);
  TEST_REAL_SIMILAR(peakShapeMetrics.start_time_at_5,2.47208309122377);
  TEST_REAL_SIMILAR(peakShapeMetrics.end_time_at_10,2.76718371294871);
  TEST_REAL_SIMILAR(peakShapeMetrics.end_time_at_5,2.75132540909723);
  TEST_REAL_SIMILAR(peakShapeMetrics.total_width,0.540983333);
  TEST_REAL_SIMILAR(peakShapeMetrics.tailing_factor,5.96347844593576);
  TEST_REAL_SIMILAR(peakShapeMetrics.asymmetry_factor,0.864852961737272);
  TEST_REAL_SIMILAR(peakShapeMetrics.baseline_delta_2_height,0.002151537878);
  TEST_REAL_SIMILAR(peakShapeMetrics.slope_of_baseline,2077);
  TEST_EQUAL(peakShapeMetrics.points_across_baseline,57);
  TEST_EQUAL(peakShapeMetrics.points_across_half_height,6);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
