// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class PeakIntegratorTest : PeakIntegrator
{
  public:
  // make protected member public
  template <typename PeakContainerConstIteratorT>
  double findPosAtPeakHeightPercent(
      PeakContainerConstIteratorT it_left,  // must not be past the end
      PeakContainerConstIteratorT it_right, // might be past the end
      PeakContainerConstIteratorT it_end,   // definitely past-the-end
      const double peak_height,
      const double percent,
      const bool is_left_half)
    {
      return findPosAtPeakHeightPercent_(it_left, it_right, it_end, peak_height, percent, is_left_half);
    }
};

START_TEST(PeakIntegrator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakIntegrator* ptr = 0;
PeakIntegrator* null_ptr = 0;

const double left = 2.472833334;
const double right = 3.022891666;

// Toy chromatogram
// data is taken from raw LC-MS/MS data points acquired for L-Glutamate in RBCs
const vector<double> position = {
  2.23095,2.239716667,2.248866667,2.25765,2.266416667,
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
  3.190133333,3.198016667,3.205916667,3.213166667
};

const vector<double> position_2 = {
  2270.93, 2272.86, 2273.16
};

const vector<double> intensity = {
  1447,2139,1699,755,1258,1070,944,1258,1573,1636,
  1762,1447,1133,1321,1762,1133,1447,2391,692,1636,2957,1321,1573,1196,1258,881,
  1384,2076,1133,1699,1384,692,1636,1133,1573,1825,1510,2391,4342,10382,17618,
  51093,153970,368094,632114,869730,962547,966489,845055,558746,417676,270942,
  184865,101619,59776,44863,31587,24036,20450,20324,11074,9879,10508,7928,7110,
  6733,6481,5726,6921,6670,5537,4971,4719,4782,5097,5789,4279,5411,4530,3524,
  2139,3335,3083,4342,4279,3083,3649,4216,4216,3964,2957,2202,2391,2643,3524,
  2328,2202,3649,2706,3020,3335,2580,2328,2894,3146,2769,2517
};

const vector<double> intensity_2 = {
  410430.0, 166125.0, 896669.0
};

const double left_past_5 = position[41];   // 2.64065
const double left_past_10 = position[42];  // 2.65125
const double left_past_50 = position[44];  // 2.672716667
const double right_past_5 = position[54];  // 2.779016667
const double right_past_10 = position[53]; // 2.768416667
const double right_past_50 = position[49]; // 2.725683333
const double left_few = position[46];      // 2.6939
const double right_few = position[48];     // 2.715083333

MSChromatogram chromatogram;
MSSpectrum spectrum;
for (Size i = 0; i < position.size(); ++i)
{
  chromatogram.push_back(ChromatogramPeak(position[i], intensity[i]));
  spectrum.push_back(Peak1D(position[i], intensity[i]));
}

MSChromatogram chromatogram_2;
MSSpectrum spectrum_2;
for (Size i = 0; i < position_2.size(); ++i)
{
  chromatogram_2.push_back(ChromatogramPeak(position_2[i], intensity_2[i]));
  spectrum_2.push_back(Peak1D(position_2[i], intensity_2[i]));
}

MSChromatogram::ConstIterator chrom_left_it = chromatogram.RTBegin(left);
MSChromatogram::ConstIterator chrom_right_it = chromatogram.RTEnd(right) - 1;
MSChromatogram::ConstIterator chrom_right_1pt_it = chromatogram.RTEnd(2.477966667) - 1;
MSChromatogram::ConstIterator chrom_right_2pt_it = chromatogram.RTEnd(2.488216667) - 1;
MSSpectrum::ConstIterator spec_left_it = spectrum.MZBegin(left);
MSSpectrum::ConstIterator spec_right_it = spectrum.MZEnd(right) - 1;
MSSpectrum::ConstIterator spec_right_1pt_it = spectrum.MZEnd(2.477966667) - 1;
MSSpectrum::ConstIterator spec_right_2pt_it = spectrum.MZEnd(2.488216667) - 1;

// To test a chromatogram with missing (5,10,50)% peak's height points
MSChromatogram::ConstIterator chrom_left_past_5_it = chromatogram.RTBegin(left_past_5);
MSChromatogram::ConstIterator chrom_right_past_5_it = chromatogram.RTEnd(right_past_5) - 1;
MSChromatogram::ConstIterator chrom_left_past_10_it = chromatogram.RTBegin(left_past_10);
MSChromatogram::ConstIterator chrom_right_past_10_it = chromatogram.RTEnd(right_past_10) - 1;
MSChromatogram::ConstIterator chrom_left_past_50_it = chromatogram.RTBegin(left_past_50);
MSChromatogram::ConstIterator chrom_right_past_50_it = chromatogram.RTEnd(right_past_50) - 1;

// To test a spectrum with missing (5,10,50)% peak's height points
MSSpectrum::ConstIterator spec_left_past_5_it = spectrum.MZBegin(left_past_5);
MSSpectrum::ConstIterator spec_right_past_5_it = spectrum.MZEnd(right_past_5) - 1;
MSSpectrum::ConstIterator spec_left_past_10_it = spectrum.MZBegin(left_past_10);
MSSpectrum::ConstIterator spec_right_past_10_it = spectrum.MZEnd(right_past_10) - 1;
MSSpectrum::ConstIterator spec_left_past_50_it = spectrum.MZBegin(left_past_50);
MSSpectrum::ConstIterator spec_right_past_50_it = spectrum.MZEnd(right_past_50) - 1;

// To test a chromatogram (and a spectrum) with few points (3 points, in this case)
MSChromatogram::ConstIterator chrom_left_few_it = chromatogram.RTBegin(left_few);
MSChromatogram::ConstIterator chrom_right_few_it = chromatogram.RTEnd(right_few) - 1;
MSSpectrum::ConstIterator spec_left_few_it = spectrum.MZBegin(left_few);
MSSpectrum::ConstIterator spec_right_few_it = spectrum.MZEnd(right_few) - 1;

constexpr const char* INTEGRATION_TYPE_INTENSITYSUM = "intensity_sum";
constexpr const char* INTEGRATION_TYPE_TRAPEZOID = "trapezoid";
constexpr const char* INTEGRATION_TYPE_SIMPSON = "simpson";
constexpr const char* BASELINE_TYPE_BASETOBASE = "base_to_base";
constexpr const char* BASELINE_TYPE_VERTICALDIVISION_MIN = "vertical_division_min";
constexpr const char* BASELINE_TYPE_VERTICALDIVISION_MAX = "vertical_division_max";


START_SECTION(PeakIntegrator())
{
  ptr = new PeakIntegrator();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PeakIntegrator())
{
  delete ptr;
}
END_SECTION

ptr = new PeakIntegrator();

START_SECTION(getParameters())
{
  Param params = ptr->getParameters();
  TEST_EQUAL(params.getValue("integration_type"), INTEGRATION_TYPE_INTENSITYSUM)
  TEST_EQUAL(params.getValue("baseline_type"), BASELINE_TYPE_BASETOBASE)
}
END_SECTION

START_SECTION(PeakBackground estimateBackground(
  const MSChromatogram& chromatogram, const double left, const double right,
  const double peak_apex_pos
) const)
{
  Param params = ptr->getParameters();
  PeakIntegrator::PeakArea pa;
  PeakIntegrator::PeakBackground pb;

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, left, right);
  pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 123446.661339019)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MIN);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, left, right);
  pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 50217)
  TEST_REAL_SIMILAR(pb.height, 881)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MAX);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, left, right);
  pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 190095)
  TEST_REAL_SIMILAR(pb.height, 3335)

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, left, right);
  pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1140.392865964)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MIN);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, left, right);
  pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 476.606316373)
  TEST_REAL_SIMILAR(pb.height, 881)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MAX);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, left, right);
  pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1804.179415555)
  TEST_REAL_SIMILAR(pb.height, 3335)
}
END_SECTION

START_SECTION(PeakBackground estimateBackground(
  const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right,
  const double peak_apex_pos
) const)
{
  Param params = ptr->getParameters();
  PeakIntegrator::PeakArea pa;
  PeakIntegrator::PeakBackground pb;

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  pb = ptr->estimateBackground(chromatogram, chrom_left_it, chrom_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 123446.661339019)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MIN);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  pb = ptr->estimateBackground(chromatogram, chrom_left_it, chrom_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 50217)
  TEST_REAL_SIMILAR(pb.height, 881)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MAX);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  pb = ptr->estimateBackground(chromatogram, chrom_left_it, chrom_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 190095)
  TEST_REAL_SIMILAR(pb.height, 3335)

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  pb = ptr->estimateBackground(chromatogram, chrom_left_it, chrom_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1140.392865964)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MIN);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  pb = ptr->estimateBackground(chromatogram, chrom_left_it, chrom_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 476.606316373)
  TEST_REAL_SIMILAR(pb.height, 881)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MAX);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  pb = ptr->estimateBackground(chromatogram, chrom_left_it, chrom_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1804.179415555)
  TEST_REAL_SIMILAR(pb.height, 3335)
}
END_SECTION

START_SECTION(PeakBackground estimateBackground(
  const MSSpectrum& spectrum, const double left, const double right,
  const double peak_apex_pos
) const)
{
  Param params = ptr->getParameters();
  PeakIntegrator::PeakArea pa;
  PeakIntegrator::PeakBackground pb;

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, left, right);
  pb = ptr->estimateBackground(spectrum, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 123446.661339019)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MIN);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, left, right);
  pb = ptr->estimateBackground(spectrum, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 50217)
  TEST_REAL_SIMILAR(pb.height, 881)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MAX);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, left, right);
  pb = ptr->estimateBackground(spectrum, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 190095)
  TEST_REAL_SIMILAR(pb.height, 3335)

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, left, right);
  pb = ptr->estimateBackground(spectrum, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1140.392865964)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MIN);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, left, right);
  pb = ptr->estimateBackground(spectrum, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 476.606316373)
  TEST_REAL_SIMILAR(pb.height, 881)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MAX);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, left, right);
  pb = ptr->estimateBackground(spectrum, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1804.179415555)
  TEST_REAL_SIMILAR(pb.height, 3335)
}
END_SECTION

START_SECTION(PeakBackground estimateBackground(
  const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right,
  const double peak_apex_pos
) const)
{
  Param params = ptr->getParameters();
  PeakIntegrator::PeakArea pa;
  PeakIntegrator::PeakBackground pb;

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  pb = ptr->estimateBackground(spectrum, spec_left_it, spec_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 123446.661339019)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MIN);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  pb = ptr->estimateBackground(spectrum, spec_left_it, spec_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 50217)
  TEST_REAL_SIMILAR(pb.height, 881)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MAX);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  pb = ptr->estimateBackground(spectrum, spec_left_it, spec_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 190095)
  TEST_REAL_SIMILAR(pb.height, 3335)

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  pb = ptr->estimateBackground(spectrum, spec_left_it, spec_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1140.392865964)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MIN);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  pb = ptr->estimateBackground(spectrum, spec_left_it, spec_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 476.606316373)
  TEST_REAL_SIMILAR(pb.height, 881)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION_MAX);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  pb = ptr->estimateBackground(spectrum, spec_left_it, spec_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1804.179415555)
  TEST_REAL_SIMILAR(pb.height, 3335)
}
END_SECTION

START_SECTION(PeakArea integratePeak(
  const MSChromatogram& chromatogram, const double left, const double right
) const)
{
  Param params = ptr->getParameters();
  PeakIntegrator::PeakArea pa;

  params.setValue("integration_type", INTEGRATION_TYPE_SIMPSON);
  ptr->setParameters(params);
  STATUS("Integration type: simpson")
  
  pa = ptr->integratePeak(chromatogram_2, 2270.93, 2273.16);
  TEST_REAL_SIMILAR(pa.area, -665788.77663627)    // TO DO: Simpson rule results in negative area for strictly positive input.
}
END_SECTION

START_SECTION(PeakArea integratePeak(
  const MSChromatogram& chromatogram, const double left, const double right
) const)
{
  Param params = ptr->getParameters();
  PeakIntegrator::PeakArea pa;
  MSChromatogram::ConstIterator it;

  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  STATUS("Integration type: intensity_sum")
  pa = ptr->integratePeak(chromatogram, left, right);
  TEST_REAL_SIMILAR(pa.area, 6768778)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  STATUS("Integration type: trapezoid")
  pa = ptr->integratePeak(chromatogram, left, right);
  TEST_REAL_SIMILAR(pa.area, 71540.2)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  STATUS("Integration type: trapezoid (1 point)")
  pa = ptr->integratePeak(chromatogram, left, 2.478);
  TEST_REAL_SIMILAR(pa.area, 0.0)
  TEST_REAL_SIMILAR(pa.height, 881.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.477966667)

  params.setValue("integration_type", INTEGRATION_TYPE_SIMPSON);
  ptr->setParameters(params);
  STATUS("Integration type: simpson (ODD number of points)")
  pa = ptr->integratePeak(chromatogram, left, right);
  TEST_REAL_SIMILAR(pa.area, 71720.443144994)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }
  STATUS("Integration type: simpson (EVEN number of points)")
  pa = ptr->integratePeak(chromatogram, left, 3.011416667); // a lower value of "right" is passed, to have 1 less point
  TEST_REAL_SIMILAR(pa.area, 71515.0792609335)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  STATUS("Integration type: simpson (1 point)")
  pa = ptr->integratePeak(chromatogram, left, 2.478);
  TEST_REAL_SIMILAR(pa.area, 0.0)
  TEST_REAL_SIMILAR(pa.height, 881.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.477966667)

  STATUS("Integration type: simpson (2 points)")
  pa = ptr->integratePeak(chromatogram, left, 2.489);
  TEST_REAL_SIMILAR(pa.area, 11.6081250000001)
  TEST_REAL_SIMILAR(pa.height, 1384)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.488216667)
}
END_SECTION

START_SECTION(PeakArea integratePeak(
  const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right
) const)
{
  Param params = ptr->getParameters();
  PeakIntegrator::PeakArea pa;
  MSChromatogram::ConstIterator it;

  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  STATUS("Integration type: intensity_sum")
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  TEST_REAL_SIMILAR(pa.area, 6768778)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  STATUS("Integration type: trapezoid")
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  TEST_REAL_SIMILAR(pa.area, 71540.2)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  STATUS("Integration type: trapezoid (1 point)")
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_1pt_it);
  TEST_REAL_SIMILAR(pa.area, 0.0)
  TEST_REAL_SIMILAR(pa.height, 881.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.477966667)

  params.setValue("integration_type", INTEGRATION_TYPE_SIMPSON);
  ptr->setParameters(params);
  STATUS("Integration type: simpson (ODD number of points)")
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  TEST_REAL_SIMILAR(pa.area, 71720.443144994)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }
  STATUS("Integration type: simpson (EVEN number of points)")
  MSChromatogram::ConstIterator chrom_right_it_less = chrom_right_it - 1;
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it_less); // a lower value of "right" is passed, to have 1 less point
  TEST_REAL_SIMILAR(pa.area, 71515.0792609335)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  STATUS("Integration type: simpson (1 point)")
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_1pt_it);
  TEST_REAL_SIMILAR(pa.area, 0.0)
  TEST_REAL_SIMILAR(pa.height, 881.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.477966667)

  STATUS("Integration type: simpson (2 points)")
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_2pt_it);
  TEST_REAL_SIMILAR(pa.area, 11.6081250000001)
  TEST_REAL_SIMILAR(pa.height, 1384)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.488216667)
}
END_SECTION

START_SECTION(PeakArea integratePeak(
  const MSSpectrum& spectrum, const double left, const double right
) const)
{
  Param params = ptr->getParameters();
  PeakIntegrator::PeakArea pa;
  MSSpectrum::ConstIterator it;

  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  STATUS("Integration type: intensity_sum")
  pa = ptr->integratePeak(spectrum, left, right);
  TEST_REAL_SIMILAR(pa.area, 6768778)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  STATUS("Integration type: trapezoid")
  pa = ptr->integratePeak(spectrum, left, right);
  TEST_REAL_SIMILAR(pa.area, 71540.2)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  STATUS("Integration type: trapezoid (1 point)")
  pa = ptr->integratePeak(spectrum, left, 2.478);
  TEST_REAL_SIMILAR(pa.area, 0.0)
  TEST_REAL_SIMILAR(pa.height, 881.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.477966667)

  params.setValue("integration_type", INTEGRATION_TYPE_SIMPSON);
  ptr->setParameters(params);
  STATUS("Integration type: simpson (ODD number of points)")
  pa = ptr->integratePeak(spectrum, left, right);
  TEST_REAL_SIMILAR(pa.area, 71720.443144994)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }
  STATUS("Integration type: simpson (EVEN number of points)")
  pa = ptr->integratePeak(spectrum, left, 3.011416667); // a lower value of "right" is passed, to have 1 less point
  TEST_REAL_SIMILAR(pa.area, 71515.0792609335)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  STATUS("Integration type: simpson (1 point)")
  pa = ptr->integratePeak(spectrum, left, 2.478);
  TEST_REAL_SIMILAR(pa.area, 0.0)
  TEST_REAL_SIMILAR(pa.height, 881.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.477966667)

  STATUS("Integration type: simpson (2 points)")
  pa = ptr->integratePeak(spectrum, left, 2.489);
  TEST_REAL_SIMILAR(pa.area, 11.6081250000001)
  TEST_REAL_SIMILAR(pa.height, 1384)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.488216667)
}
END_SECTION

START_SECTION(PeakArea integratePeak(
  const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right
) const)
{
  Param params = ptr->getParameters();
  PeakIntegrator::PeakArea pa;
  MSSpectrum::ConstIterator it;

  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  STATUS("Integration type: intensity_sum")
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  TEST_REAL_SIMILAR(pa.area, 6768778)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  STATUS("Integration type: trapezoid")
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  TEST_REAL_SIMILAR(pa.area, 71540.2)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  STATUS("Integration type: trapezoid (1 point)")
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_1pt_it);
  TEST_REAL_SIMILAR(pa.area, 0.0)
  TEST_REAL_SIMILAR(pa.height, 881.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.477966667)

  params.setValue("integration_type", INTEGRATION_TYPE_SIMPSON);
  ptr->setParameters(params);
  STATUS("Integration type: simpson (ODD number of points)")
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  TEST_REAL_SIMILAR(pa.area, 71720.443144994)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }
  STATUS("Integration type: simpson (EVEN number of points)")
  MSSpectrum::ConstIterator spec_right_it_less = spec_right_it - 1;
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it_less); // a lower value of "right" is passed, to have 1 less point
  TEST_REAL_SIMILAR(pa.area, 71515.0792609335)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); i += 4, it +=4)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  STATUS("Integration type: simpson (1 point)")
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_1pt_it);
  TEST_REAL_SIMILAR(pa.area, 0.0)
  TEST_REAL_SIMILAR(pa.height, 881.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.477966667)

  STATUS("Integration type: simpson (2 points)")
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_2pt_it);
  TEST_REAL_SIMILAR(pa.area, 11.6081250000001)
  TEST_REAL_SIMILAR(pa.height, 1384)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.488216667)
}
END_SECTION

START_SECTION(PeakShapeMetrics calculatePeakShapeMetrics(
  const MSChromatogram& chromatogram, const double left, const double right,
  const double peak_height, const double peak_apex_pos
) const)
{
  PeakIntegrator::PeakArea pa;
  PeakIntegrator::PeakShapeMetrics psm;
  pa = ptr->integratePeak(chromatogram, left, right);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, left, right, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.width_at_5, 0.15955)
  TEST_REAL_SIMILAR(psm.width_at_10, 0.138366667)
  TEST_REAL_SIMILAR(psm.width_at_50, 0.0741500000000004)
  TEST_REAL_SIMILAR(psm.start_position_at_5, 2.630066667)
  TEST_REAL_SIMILAR(psm.start_position_at_10, 2.64065)
  TEST_REAL_SIMILAR(psm.start_position_at_50, 2.662116667)
  TEST_REAL_SIMILAR(psm.end_position_at_5, 2.789616667)
  TEST_REAL_SIMILAR(psm.end_position_at_10, 2.779016667)
  TEST_REAL_SIMILAR(psm.end_position_at_50, 2.736266667)
  TEST_REAL_SIMILAR(psm.total_width, 0.540983333)
  TEST_REAL_SIMILAR(psm.tailing_factor, 1.07176444725376)
  TEST_REAL_SIMILAR(psm.asymmetry_factor, 1.16705821456539)
  TEST_REAL_SIMILAR(psm.slope_of_baseline, 2454)
  TEST_REAL_SIMILAR(psm.baseline_delta_2_height, 0.00253908735640033)
  TEST_EQUAL(psm.points_across_baseline, 57)
  TEST_EQUAL(psm.points_across_half_height, 6)
  pa = ptr->integratePeak(chromatogram, left_past_5, right_past_5);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, left_past_5, right_past_5, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_5)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_5)
  pa = ptr->integratePeak(chromatogram, left_past_10, right_past_10);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, left_past_10, right_past_10, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_10)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_10)
  TEST_REAL_SIMILAR(psm.start_position_at_10, left_past_10)
  TEST_REAL_SIMILAR(psm.end_position_at_10, right_past_10)
  pa = ptr->integratePeak(chromatogram, left_past_50, right_past_50);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, left_past_50, right_past_50, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_50)
  TEST_REAL_SIMILAR(psm.start_position_at_10, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_10, right_past_50)
  TEST_REAL_SIMILAR(psm.start_position_at_50, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_50, right_past_50)
  pa = ptr->integratePeak(chromatogram, left_few, right_few);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, left_few, right_few, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_few)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_few)
}
END_SECTION

START_SECTION(PeakShapeMetrics calculatePeakShapeMetrics(
  const MSChromatogram& chromatogram, MSChromatogram::ConstIterator& left, MSChromatogram::ConstIterator& right,
  const double peak_height, const double peak_apex_pos
) const)
{
  PeakIntegrator::PeakArea pa;
  PeakIntegrator::PeakShapeMetrics psm;
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, chrom_left_it, chrom_right_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.width_at_5, 0.15955)
  TEST_REAL_SIMILAR(psm.width_at_10, 0.138366667)
  TEST_REAL_SIMILAR(psm.width_at_50, 0.0741500000000004)
  TEST_REAL_SIMILAR(psm.start_position_at_5, 2.630066667)
  TEST_REAL_SIMILAR(psm.start_position_at_10, 2.64065)
  TEST_REAL_SIMILAR(psm.start_position_at_50, 2.662116667)
  TEST_REAL_SIMILAR(psm.end_position_at_5, 2.789616667)
  TEST_REAL_SIMILAR(psm.end_position_at_10, 2.779016667)
  TEST_REAL_SIMILAR(psm.end_position_at_50, 2.736266667)
  TEST_REAL_SIMILAR(psm.total_width, 0.540983333)
  TEST_REAL_SIMILAR(psm.tailing_factor, 1.07176444725376)
  TEST_REAL_SIMILAR(psm.asymmetry_factor, 1.16705821456539)
  TEST_REAL_SIMILAR(psm.slope_of_baseline, 2454)
  TEST_REAL_SIMILAR(psm.baseline_delta_2_height, 0.00253908735640033)
  TEST_EQUAL(psm.points_across_baseline, 57)
  TEST_EQUAL(psm.points_across_half_height, 6)
  pa = ptr->integratePeak(chromatogram, chrom_left_past_5_it, chrom_right_past_5_it);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, chrom_left_past_5_it, chrom_right_past_5_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_5)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_5)
  pa = ptr->integratePeak(chromatogram, chrom_left_past_10_it, chrom_right_past_10_it);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, chrom_left_past_10_it, chrom_right_past_10_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_10)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_10)
  TEST_REAL_SIMILAR(psm.start_position_at_10, left_past_10)
  TEST_REAL_SIMILAR(psm.end_position_at_10, right_past_10)
  pa = ptr->integratePeak(chromatogram, chrom_left_past_50_it, chrom_right_past_50_it);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, chrom_left_past_50_it, chrom_right_past_50_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_50)
  TEST_REAL_SIMILAR(psm.start_position_at_10, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_10, right_past_50)
  TEST_REAL_SIMILAR(psm.start_position_at_50, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_50, right_past_50)
  pa = ptr->integratePeak(chromatogram, chrom_left_few_it, chrom_right_few_it);
  psm = ptr->calculatePeakShapeMetrics(chromatogram, chrom_left_few_it, chrom_right_few_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_few)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_few)
}
END_SECTION

START_SECTION(PeakShapeMetrics calculatePeakShapeMetrics(
  const MSSpectrum& spectrum, const double left, const double right,
  const double peak_height, const double peak_apex_pos
) const)
{
  PeakIntegrator::PeakArea pa;
  PeakIntegrator::PeakShapeMetrics psm;
  pa = ptr->integratePeak(spectrum, left, right);
  psm = ptr->calculatePeakShapeMetrics(spectrum, left, right, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.width_at_5, 0.15955)
  TEST_REAL_SIMILAR(psm.width_at_10, 0.138366667)
  TEST_REAL_SIMILAR(psm.width_at_50, 0.0741500000000004)
  TEST_REAL_SIMILAR(psm.start_position_at_5, 2.630066667)
  TEST_REAL_SIMILAR(psm.start_position_at_10, 2.64065)
  TEST_REAL_SIMILAR(psm.start_position_at_50, 2.662116667)
  TEST_REAL_SIMILAR(psm.end_position_at_5, 2.789616667)
  TEST_REAL_SIMILAR(psm.end_position_at_10, 2.779016667)
  TEST_REAL_SIMILAR(psm.end_position_at_50, 2.736266667)
  TEST_REAL_SIMILAR(psm.total_width, 0.540983333)
  TEST_REAL_SIMILAR(psm.tailing_factor, 1.07176444725376)
  TEST_REAL_SIMILAR(psm.asymmetry_factor, 1.16705821456539)
  TEST_REAL_SIMILAR(psm.slope_of_baseline, 2454)
  TEST_REAL_SIMILAR(psm.baseline_delta_2_height, 0.00253908735640033)
  TEST_EQUAL(psm.points_across_baseline, 57)
  TEST_EQUAL(psm.points_across_half_height, 6)
  pa = ptr->integratePeak(spectrum, left_past_5, right_past_5);
  psm = ptr->calculatePeakShapeMetrics(spectrum, left_past_5, right_past_5, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_5)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_5)
  pa = ptr->integratePeak(spectrum, left_past_10, right_past_10);
  psm = ptr->calculatePeakShapeMetrics(spectrum, left_past_10, right_past_10, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_10)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_10)
  TEST_REAL_SIMILAR(psm.start_position_at_10, left_past_10)
  TEST_REAL_SIMILAR(psm.end_position_at_10, right_past_10)
  pa = ptr->integratePeak(spectrum, left_past_50, right_past_50);
  psm = ptr->calculatePeakShapeMetrics(spectrum, left_past_50, right_past_50, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_50)
  TEST_REAL_SIMILAR(psm.start_position_at_10, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_10, right_past_50)
  TEST_REAL_SIMILAR(psm.start_position_at_50, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_50, right_past_50)
  pa = ptr->integratePeak(spectrum, left_few, right_few);
  psm = ptr->calculatePeakShapeMetrics(spectrum, left_few, right_few, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_few)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_few)
}
END_SECTION

START_SECTION(PeakShapeMetrics calculatePeakShapeMetrics(
  const MSSpectrum& spectrum, MSSpectrum::ConstIterator& left, MSSpectrum::ConstIterator& right,
  const double peak_height, const double peak_apex_pos
) const)
{
  PeakIntegrator::PeakArea pa;
  PeakIntegrator::PeakShapeMetrics psm;
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  psm = ptr->calculatePeakShapeMetrics(spectrum, spec_left_it, spec_right_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.width_at_5, 0.15955)
  TEST_REAL_SIMILAR(psm.width_at_10, 0.138366667)
  TEST_REAL_SIMILAR(psm.width_at_50, 0.0741500000000004)
  TEST_REAL_SIMILAR(psm.start_position_at_5, 2.630066667)
  TEST_REAL_SIMILAR(psm.start_position_at_10, 2.64065)
  TEST_REAL_SIMILAR(psm.start_position_at_50, 2.662116667)
  TEST_REAL_SIMILAR(psm.end_position_at_5, 2.789616667)
  TEST_REAL_SIMILAR(psm.end_position_at_10, 2.779016667)
  TEST_REAL_SIMILAR(psm.end_position_at_50, 2.736266667)
  TEST_REAL_SIMILAR(psm.total_width, 0.540983333)
  TEST_REAL_SIMILAR(psm.tailing_factor, 1.07176444725376)
  TEST_REAL_SIMILAR(psm.asymmetry_factor, 1.16705821456539)
  TEST_REAL_SIMILAR(psm.slope_of_baseline, 2454)
  TEST_REAL_SIMILAR(psm.baseline_delta_2_height, 0.00253908735640033)
  TEST_EQUAL(psm.points_across_baseline, 57)
  TEST_EQUAL(psm.points_across_half_height, 6)
  pa = ptr->integratePeak(spectrum, spec_left_past_5_it, spec_right_past_5_it);
  psm = ptr->calculatePeakShapeMetrics(spectrum, spec_left_past_5_it, spec_right_past_5_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_5)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_5)
  pa = ptr->integratePeak(spectrum, spec_left_past_10_it, spec_right_past_10_it);
  psm = ptr->calculatePeakShapeMetrics(spectrum, spec_left_past_10_it, spec_right_past_10_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_10)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_10)
  TEST_REAL_SIMILAR(psm.start_position_at_10, left_past_10)
  TEST_REAL_SIMILAR(psm.end_position_at_10, right_past_10)
  pa = ptr->integratePeak(spectrum, spec_left_past_50_it, spec_right_past_50_it);
  psm = ptr->calculatePeakShapeMetrics(spectrum, spec_left_past_50_it, spec_right_past_50_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_past_50)
  TEST_REAL_SIMILAR(psm.start_position_at_10, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_10, right_past_50)
  TEST_REAL_SIMILAR(psm.start_position_at_50, left_past_50)
  TEST_REAL_SIMILAR(psm.end_position_at_50, right_past_50)
  pa = ptr->integratePeak(spectrum, spec_left_few_it, spec_right_few_it);
  psm = ptr->calculatePeakShapeMetrics(spectrum, spec_left_few_it, spec_right_few_it, pa.height, pa.apex_pos);
  TEST_REAL_SIMILAR(psm.start_position_at_5, left_few)
  TEST_REAL_SIMILAR(psm.end_position_at_5, right_few)
}
END_SECTION

START_SECTION([EXTRA]  template <typename PeakContainerConstIteratorT> double findPosAtPeakHeightPercent_(...))
{
  PeakIntegratorTest pit;
  double pos = pit.findPosAtPeakHeightPercent(spectrum.begin(), spectrum.end(), spectrum.end(), 0.0, 0.0, true);
  TEST_EQUAL(pos, spectrum[0].getPos()); // find first non-zero peak

  pos = pit.findPosAtPeakHeightPercent(spectrum.begin(), spectrum.end(), spectrum.end(), 0.0, 0.0, false);
  TEST_EQUAL(pos, spectrum.back().getPos()); // find non-zero peak from end

  // corner cases: just a single point in range
  pos = pit.findPosAtPeakHeightPercent(spectrum.begin(), spectrum.begin() + 1, spectrum.begin() + 1, 0.0, 0.0, false);
  TEST_EQUAL(pos, spectrum[0].getPos()); // return the only peak there is
  pos = pit.findPosAtPeakHeightPercent(spectrum.begin(), spectrum.begin() + 1, spectrum.begin() + 1, 0.0, 0.0, true);
  TEST_EQUAL(pos, spectrum[0].getPos()); // return the only peak there is


  // corner cases: empty range
  TEST_EXCEPTION(Exception::InvalidRange, pit.findPosAtPeakHeightPercent(spectrum.end(), spectrum.end(), spectrum.end(), 0.0, 0.0, false))
  TEST_EXCEPTION(Exception::InvalidRange, pit.findPosAtPeakHeightPercent(spectrum.end(), spectrum.end(), spectrum.end(), 0.0, 0.0, true))
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
