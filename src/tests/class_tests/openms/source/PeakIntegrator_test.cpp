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

START_TEST(PeakIntegrator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakIntegrator* ptr = 0;
PeakIntegrator* null_ptr = 0;

const double left = 2.472833334;
const double right = 3.022891666;

// Toy chromatogram
// data is taken from raw LC-MS/MS data points acquired for L-Glutamate in RBCs
vector<double> position = {
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

vector<double> intensity = {
  1447,2139,1699,755,1258,1070,944,1258,1573,1636,
  1762,1447,1133,1321,1762,1133,1447,2391,692,1636,2957,1321,1573,1196,1258,881,
  1384,2076,1133,1699,1384,692,1636,1133,1573,1825,1510,2391,4342,10382,17618,
  51093,153970,368094,632114,869730,962547,966489,845055,558746,417676,270942,
  184865,101619,59776,44863,31587,24036,20450,20324,11074,9879,10508,7928,7110,
  6733,6481,5726,6921,6670,5537,4971,4719,4782,5097,5789,4279,5411,4530,3524,
  2139,3335,3083,4342,4279,3083,3649,4216,4216,3964,2957,2202,2391,2643,3524,
  2328,2202,3649,2706,3020,3335,2580,2328,2894,3146,2769,2517
};

MSChromatogram chromatogram;
MSSpectrum spectrum;
for (Size i = 0; i < position.size(); ++i)
{
  chromatogram.push_back(ChromatogramPeak(position[i], intensity[i]));
  spectrum.push_back(Peak1D(position[i], intensity[i]));
}

MSChromatogram::ConstIterator chrom_left_it = chromatogram.RTBegin(left);
MSChromatogram::ConstIterator chrom_right_it = chromatogram.RTEnd(right) - 1;
MSSpectrum::ConstIterator spec_left_it = spectrum.MZBegin(left);
MSSpectrum::ConstIterator spec_right_it = spectrum.MZEnd(right) - 1;

constexpr const char* INTEGRATION_TYPE_INTENSITYSUM = "intensity_sum";
constexpr const char* INTEGRATION_TYPE_TRAPEZOID = "trapezoid";
constexpr const char* INTEGRATION_TYPE_SIMPSON = "simpson";
constexpr const char* BASELINE_TYPE_BASETOBASE = "base_to_base";
constexpr const char* BASELINE_TYPE_VERTICALDIVISION = "vertical_division";

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

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, left, right);
  pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 50217)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, left, right);
  pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1140.392865964)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, left, right);
  pb = ptr->estimateBackground(chromatogram, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 476.606316373)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)
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

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  pb = ptr->estimateBackground(chromatogram, chrom_left_it, chrom_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 50217)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  pb = ptr->estimateBackground(chromatogram, chrom_left_it, chrom_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1140.392865964)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  pb = ptr->estimateBackground(chromatogram, chrom_left_it, chrom_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 476.606316373)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)
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

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, left, right);
  pb = ptr->estimateBackground(spectrum, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 50217)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, left, right);
  pb = ptr->estimateBackground(spectrum, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1140.392865964)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, left, right);
  pb = ptr->estimateBackground(spectrum, left, right, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 476.606316373)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)
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

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION);
  params.setValue("integration_type", INTEGRATION_TYPE_INTENSITYSUM);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  pb = ptr->estimateBackground(spectrum, spec_left_it, spec_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 50217)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_BASETOBASE);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  pb = ptr->estimateBackground(spectrum, spec_left_it, spec_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 1140.392865964)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)

  params.setValue("baseline_type", BASELINE_TYPE_VERTICALDIVISION);
  params.setValue("integration_type", INTEGRATION_TYPE_TRAPEZOID);
  ptr->setParameters(params);
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  pb = ptr->estimateBackground(spectrum, spec_left_it, spec_right_it, pa.apex_pos);
  TEST_REAL_SIMILAR(pb.area, 476.606316373)
  TEST_REAL_SIMILAR(pb.height, 1908.59690598823)
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  params.setValue("integration_type", INTEGRATION_TYPE_SIMPSON);
  ptr->setParameters(params);
  STATUS("Integration type: simpson (ODD number of points)")
  pa = ptr->integratePeak(chromatogram, left, right);
  TEST_REAL_SIMILAR(pa.area, 71720.443144994)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  params.setValue("integration_type", INTEGRATION_TYPE_SIMPSON);
  ptr->setParameters(params);
  STATUS("Integration type: simpson (ODD number of points)")
  pa = ptr->integratePeak(chromatogram, chrom_left_it, chrom_right_it);
  TEST_REAL_SIMILAR(pa.area, 71720.443144994)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = chromatogram.RTBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getRT())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  params.setValue("integration_type", INTEGRATION_TYPE_SIMPSON);
  ptr->setParameters(params);
  STATUS("Integration type: simpson (ODD number of points)")
  pa = ptr->integratePeak(spectrum, left, right);
  TEST_REAL_SIMILAR(pa.area, 71720.443144994)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }

  params.setValue("integration_type", INTEGRATION_TYPE_SIMPSON);
  ptr->setParameters(params);
  STATUS("Integration type: simpson (ODD number of points)")
  pa = ptr->integratePeak(spectrum, spec_left_it, spec_right_it);
  TEST_REAL_SIMILAR(pa.area, 71720.443144994)
  TEST_REAL_SIMILAR(pa.height, 966489.0)
  TEST_REAL_SIMILAR(pa.apex_pos, 2.7045)
  it = spectrum.MZBegin(left);
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
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
  for (Size i = 0; i < pa.hull_points.size(); ++i, ++it)
  {
    TEST_REAL_SIMILAR(pa.hull_points[i][0], it->getMZ())
    TEST_REAL_SIMILAR(pa.hull_points[i][1], it->getIntensity())
  }
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
  TEST_REAL_SIMILAR(psm.width_at_5, 0.231263425125414)
  TEST_REAL_SIMILAR(psm.width_at_10, 0.134762234301732)
  TEST_REAL_SIMILAR(psm.width_at_50, 0.0595791540757924)
  TEST_REAL_SIMILAR(psm.start_position_at_5, 2.51268515480125)
  TEST_REAL_SIMILAR(psm.start_position_at_10, 2.63222565817823)
  TEST_REAL_SIMILAR(psm.start_position_at_50, 2.65391757114759)
  TEST_REAL_SIMILAR(psm.end_position_at_5, 2.74394857992666)
  TEST_REAL_SIMILAR(psm.end_position_at_10, 2.76698789247996)
  TEST_REAL_SIMILAR(psm.end_position_at_50, 2.71349672522338)
  TEST_REAL_SIMILAR(psm.total_width, 0.540983333)
  TEST_REAL_SIMILAR(psm.tailing_factor, 5.86240177860251)
  TEST_REAL_SIMILAR(psm.asymmetry_factor, 0.864593034054243)
  TEST_REAL_SIMILAR(psm.slope_of_baseline, 2454)
  TEST_REAL_SIMILAR(psm.baseline_delta_2_height, 0.00253908735640033)
  TEST_EQUAL(psm.points_across_baseline, 57)
  TEST_EQUAL(psm.points_across_half_height, 6)
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
  TEST_REAL_SIMILAR(psm.width_at_5, 0.231263425125414)
  TEST_REAL_SIMILAR(psm.width_at_10, 0.134762234301732)
  TEST_REAL_SIMILAR(psm.width_at_50, 0.0595791540757924)
  TEST_REAL_SIMILAR(psm.start_position_at_5, 2.51268515480125)
  TEST_REAL_SIMILAR(psm.start_position_at_10, 2.63222565817823)
  TEST_REAL_SIMILAR(psm.start_position_at_50, 2.65391757114759)
  TEST_REAL_SIMILAR(psm.end_position_at_5, 2.74394857992666)
  TEST_REAL_SIMILAR(psm.end_position_at_10, 2.76698789247996)
  TEST_REAL_SIMILAR(psm.end_position_at_50, 2.71349672522338)
  TEST_REAL_SIMILAR(psm.total_width, 0.540983333)
  TEST_REAL_SIMILAR(psm.tailing_factor, 5.86240177860251)
  TEST_REAL_SIMILAR(psm.asymmetry_factor, 0.864593034054243)
  TEST_REAL_SIMILAR(psm.slope_of_baseline, 2454)
  TEST_REAL_SIMILAR(psm.baseline_delta_2_height, 0.00253908735640033)
  TEST_EQUAL(psm.points_across_baseline, 57)
  TEST_EQUAL(psm.points_across_half_height, 6)
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
  TEST_REAL_SIMILAR(psm.width_at_5, 0.231263425125414)
  TEST_REAL_SIMILAR(psm.width_at_10, 0.134762234301732)
  TEST_REAL_SIMILAR(psm.width_at_50, 0.0595791540757924)
  TEST_REAL_SIMILAR(psm.start_position_at_5, 2.51268515480125)
  TEST_REAL_SIMILAR(psm.start_position_at_10, 2.63222565817823)
  TEST_REAL_SIMILAR(psm.start_position_at_50, 2.65391757114759)
  TEST_REAL_SIMILAR(psm.end_position_at_5, 2.74394857992666)
  TEST_REAL_SIMILAR(psm.end_position_at_10, 2.76698789247996)
  TEST_REAL_SIMILAR(psm.end_position_at_50, 2.71349672522338)
  TEST_REAL_SIMILAR(psm.total_width, 0.540983333)
  TEST_REAL_SIMILAR(psm.tailing_factor, 5.86240177860251)
  TEST_REAL_SIMILAR(psm.asymmetry_factor, 0.864593034054243)
  TEST_REAL_SIMILAR(psm.slope_of_baseline, 2454)
  TEST_REAL_SIMILAR(psm.baseline_delta_2_height, 0.00253908735640033)
  TEST_EQUAL(psm.points_across_baseline, 57)
  TEST_EQUAL(psm.points_across_half_height, 6)
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
  TEST_REAL_SIMILAR(psm.width_at_5, 0.231263425125414)
  TEST_REAL_SIMILAR(psm.width_at_10, 0.134762234301732)
  TEST_REAL_SIMILAR(psm.width_at_50, 0.0595791540757924)
  TEST_REAL_SIMILAR(psm.start_position_at_5, 2.51268515480125)
  TEST_REAL_SIMILAR(psm.start_position_at_10, 2.63222565817823)
  TEST_REAL_SIMILAR(psm.start_position_at_50, 2.65391757114759)
  TEST_REAL_SIMILAR(psm.end_position_at_5, 2.74394857992666)
  TEST_REAL_SIMILAR(psm.end_position_at_10, 2.76698789247996)
  TEST_REAL_SIMILAR(psm.end_position_at_50, 2.71349672522338)
  TEST_REAL_SIMILAR(psm.total_width, 0.540983333)
  TEST_REAL_SIMILAR(psm.tailing_factor, 5.86240177860251)
  TEST_REAL_SIMILAR(psm.asymmetry_factor, 0.864593034054243)
  TEST_REAL_SIMILAR(psm.slope_of_baseline, 2454)
  TEST_REAL_SIMILAR(psm.baseline_delta_2_height, 0.00253908735640033)
  TEST_EQUAL(psm.points_across_baseline, 57)
  TEST_EQUAL(psm.points_across_half_height, 6)
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
