// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>

#include <boost/assign/std/vector.hpp>

#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

using namespace OpenMS;
using namespace std;

typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;
// TODO also test the picker with the LightTransition interface
// typedef MRMTransitionGroup<RichPeakChromatogram, OpenSwath::LightTransition> MRMTransitionGroupType;
typedef MRMTransitionGroup<RichPeakChromatogram, ReactionMonitoringTransition> MRMTransitionGroupType;

void setup_transition_group(MRMTransitionGroupType & transition_group)
{
  // this is a simulated SRM experiment where the two traces are not sampled at
  // the exact same time points, thus a resampling is necessary before applying
  // the algorithm.
  static const double rtdata_1[] = {1474.34, 1477.11, 1479.88, 1482.64, 1485.41, 1488.19, 1490.95, 1493.72, 1496.48, 1499.25, 1502.03, 1504.8 , 1507.56, 1510.33, 1513.09, 1515.87, 1518.64, 1521.42};
  static const double rtdata_2[] = {1473.55, 1476.31, 1479.08, 1481.84, 1484.61, 1487.39, 1490.15, 1492.92, 1495.69, 1498.45, 1501.23, 1504   , 1506.76, 1509.53, 1512.29, 1515.07, 1517.84, 1520.62};

  static const double intdata_1[] = {3.26958, 3.74189, 3.31075, 86.1901, 3.47528, 387.864, 13281  , 6375.84, 39852.6, 2.66726, 612.747, 3.34313, 793.12 , 3.29156, 4.00586, 4.1591 , 3.23035, 3.90591};
  static const double intdata_2[] =  {3.44054 , 2142.31 , 3.58763 , 3076.97 , 6663.55 , 45681   , 157694  , 122844  , 86034.7 , 85391.1 , 15992.8 , 2293.94 , 6934.85 , 2735.18 , 459.413 , 3.93863 , 3.36564 , 3.44005};

  {
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
  }

  {
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
  }
}

START_TEST(MRMTransitionGroupPicker, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMTransitionGroupPicker* ptr = 0;
MRMTransitionGroupPicker* nullPointer = 0;

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
  ConvexHull2D::PointArrayType data1_points = mrmfeature.getFeature("2").getConvexHulls()[0].getHullPoints();
  double sum = 0.0;
  for (ConvexHull2D::PointArrayType::iterator it = data1_points.begin(); it != data1_points.end(); it++)
  {
    sum += it->getY();
  }
  TEST_REAL_SIMILAR(sum, 507385.32);
  TEST_REAL_SIMILAR(mrmfeature.getFeature("2").getIntensity(), 507385.32);
  ConvexHull2D::PointArrayType data2_points = mrmfeature.getFeature("1").getConvexHulls()[0].getHullPoints();
  sum = 0.0;
  for (ConvexHull2D::PointArrayType::iterator it = data2_points.begin(); it != data2_points.end(); it++)
  {
    sum += it->getY();
  }
  TEST_REAL_SIMILAR(sum, 59989.8287208466);
  TEST_REAL_SIMILAR(mrmfeature.getFeature("1").getIntensity(), 59989.8287208466);
}
END_SECTION

START_SECTION(void pickChromatogram(const RichPeakChromatogram &chromatogram, RichPeakChromatogram &smoothed_chrom, RichPeakChromatogram &picked_chrom))
{
  MRMTransitionGroupType transition_group;
  setup_transition_group(transition_group);
  std::vector< RichPeakChromatogram > picked_chroms;

  RichPeakChromatogram picked_chrom, smoothed_chrom;
  RichPeakChromatogram chrom = transition_group.getChromatograms()[0];
  MRMTransitionGroupPicker picker;
  picker.pickChromatogram(chrom, smoothed_chrom, picked_chrom);

  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), 3);

  // Peak picking is done on the smoothed data by cubic spline interpolation
  // and searching for the point with zero derivative.
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 9981.76460102146);
  TEST_REAL_SIMILAR( picked_chrom[0].getMZ(), 1495.11321013749);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[0][0], 59509.4); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[1][0], 1490.95); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[2][0], 1496.48); // rightWidth

  chrom = transition_group.getChromatograms()[1];
  picker.pickChromatogram(chrom, smoothed_chrom, picked_chrom);

  TEST_EQUAL( picked_chrom.size(), 1);
  TEST_EQUAL( picked_chrom.getFloatDataArrays().size(), 3);

  // Peak picking is done on the smoothed data by cubic spline interpolation
  // and searching for the point with zero derivative.
  TEST_REAL_SIMILAR( picked_chrom[0].getIntensity(), 78719.134569503);
  TEST_REAL_SIMILAR( picked_chrom[0].getMZ(), 1492.830608593);
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[0][0], 523378); // IntegratedIntensity
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[1][0], 1481.84); // leftWidth
  TEST_REAL_SIMILAR( picked_chrom.getFloatDataArrays()[2][0], 1501.23); // rightWidth
}
END_SECTION

START_SECTION( ( template < typename SpectrumT, typename TransitionT > MRMFeature createMRMFeature(MRMTransitionGroup< SpectrumT, TransitionT > &transition_group, std::vector< SpectrumT > &picked_chroms, int &chr_idx, int &peak_idx)))
{
  MRMTransitionGroupType transition_group;
  setup_transition_group(transition_group);
  std::vector< RichPeakChromatogram > picked_chroms;

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

    picked_chroms.push_back(picked_chrom);
  }

  // create the corresponding first mrm feature
  int chr_idx = 1, peak_idx = 0;
  MRMTransitionGroupPicker picker;
  MRMFeature mrmfeature = picker.createMRMFeature(transition_group, picked_chroms, chr_idx, peak_idx);
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
  // peaks in the first chromatogram and and the first peak in the second
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

