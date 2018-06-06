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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerMaxima.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <cmath>
#define PI 3.141592653589793

using namespace OpenMS;
using namespace std;


std::vector<PeakPickerMaxima::PeakCandidate> ppmax_pick(MSSpectrum& spec, PeakPickerMaxima& pp_max)
{
  std::vector<PeakPickerMaxima::PeakCandidate> pc;
  std::vector<double> mz_array(spec.size()), int_array(spec.size());
  for (Size p = 0; p < spec.size(); ++p)
  {
    mz_array[p] = spec[p].getMZ();
    int_array[p] = spec[p].getIntensity();
  }
  pp_max.pick(mz_array, int_array, pc);
  return pc;
}

double getGauss(double mu, double sigma, double x)
{
  return (1.0/(sigma*sqrt(2*PI)) * exp( -(x-mu)*(x-mu)/(2.0*sigma*sigma)));
}

void generateTestData(std::vector<double>& x, std::vector<double>& y, double deltax, double int_multiplicator)
{
  for (Size i = 0; i < 20; i++)
  {
    x.push_back(i + deltax);
    y.push_back( getGauss(10,5,i)*int_multiplicator );
  }
}

START_TEST(PeakPickerMaxima, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerMaxima* ptr = nullptr;
PeakPickerMaxima* nullPointer = nullptr;
START_SECTION((PeakPickerMaxima()))
  ptr = new PeakPickerMaxima(0,0,0);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeakPickerMaxima()))
  delete ptr;
END_SECTION


  /*
   * Python code:
   *

import math, numpy
def gauss(mu, sigma, x):
  return (1.0/(sigma*math.sqrt(2*numpy.pi)) * math.exp( -(x- mu)**2.0/(2.0*sigma*sigma)))

y = [ gauss(10.0,5.0,i) for i in range(20) ]
x = [ i for i in range(20) ]
xx = [x[2*i+1] for i in range(10)]
yy = [y[2*i+1] for i in range(10)]
xrand = [xit + 0.5*(random.random()-0.5) for xit in x]
yrand = [yit + max(y)*0.25*(random.random()-0.5) for yit in y]


  max(y)
  0.07978845608028655

  y[9]
  0.07820853879509118
  */

static const double arr1[] = {
  -0.04732030993393693,
  0.8914924331847927,
  2.242251028116535,
  2.8489997501981357,
  4.1663063904956,
  4.770450183181009,
  5.8026362378461815,
  7.067111623628946,
  7.968908908421478,
  8.959060860876802,
  10.005216076641757,
  11.196610814166815,
  12.116982813029852,
  12.977162356375791,
  14.100893620425213,
  15.222529820550236,
  15.890138455823378,
  16.771297077874447,
  18.074575078568163,
  19.093826607410218 
};
std::vector<double> xrand (arr1, arr1 + sizeof(arr1) / sizeof(arr1[0]) );
static const double arr2[] = {
  0.007050921402452849,
  0.006089927970860897,
  0.016432781047452296,
  0.0351300895434513,
  0.03593042409081977,
  0.0415877855954923,
  0.054303625272399056,
  0.061883694314788226,
  0.07837041224348473,
  0.07652346739035985,
  0.07886053568902987,
  0.07959292444993592,
  0.0708551147646475,
  0.05812243127463338,
  0.06262580825607922,
  0.046054457061387874,
  0.046166241756351346,
  0.023122371466895074,
  0.03182678605750754,
  0.009819277289325083

};
std::vector<double> yrand (arr2, arr2 + sizeof(arr2) / sizeof(arr2[0]) );

START_SECTION([EXTRA](pick single peak))
{
  PeakPickerMaxima pp_max(0.0);
  std::vector<double> mz_array, int_array;
  generateTestData(mz_array, int_array, 0, 1.0);

  // Test Gaussian function
  {
  std::vector<PeakPickerMaxima::PeakCandidate> pc;
  pp_max.pick(mz_array, int_array, pc);

  TEST_EQUAL(pc.size(), 1)
  TEST_EQUAL(pc[0].pos, 10)
  TEST_REAL_SIMILAR(pc[0].int_max, 0.07978845608028655)
  TEST_REAL_SIMILAR(pc[0].mz_max, 10)

  }

  // Test Gaussian function with a zero left
  {
  std::vector<PeakPickerMaxima::PeakCandidate> pc;
  std::vector<double> mz_array_n(mz_array), int_array_n(int_array);
  int_array_n[8] = 0.0;
  pp_max.pick(mz_array_n, int_array_n, pc);

  TEST_EQUAL(pc.size(), 1)
  TEST_EQUAL(pc[0].pos, 10)
  TEST_EQUAL(pc[0].left_boundary, 8)
  TEST_EQUAL(pc[0].right_boundary, 19)

  }

  // Test Gaussian function with a zero right
  {
  std::vector<PeakPickerMaxima::PeakCandidate> pc;
  std::vector<double> mz_array_n(mz_array), int_array_n(int_array);
  int_array_n[12] = 0.0;
  pp_max.pick(mz_array_n, int_array_n, pc);

  TEST_EQUAL(pc.size(), 1)
  TEST_EQUAL(pc[0].pos, 10)
  TEST_EQUAL(pc[0].left_boundary, 0)
  TEST_EQUAL(pc[0].right_boundary, 12)

  }

  // Re-sample at every second point
  {
    std::vector<double> mz_array_mut, int_array_mut;
    for (Size i = 0; i < 10; i++)
    {
      mz_array_mut.push_back( mz_array[2*i+1]);
      int_array_mut.push_back( int_array[2*i+1]);
    }

    std::vector<PeakPickerMaxima::PeakCandidate> pc;
    int_array_mut[4] += 0.0001; // needs a small delta to be backwards compatible
    pp_max.pick(mz_array_mut, int_array_mut, pc);

    TEST_EQUAL(pc.size(), 1)
    TEST_EQUAL(pc[0].pos, 4)
    TOLERANCE_RELATIVE(1.005);
    TEST_REAL_SIMILAR(pc[0].int_max, 0.07978845608028655)
    TEST_REAL_SIMILAR(pc[0].mz_max, 10)
  }

  // Re-sample at every second point
  // Introduce a small m/z error
  {
    std::vector<PeakPickerMaxima::PeakCandidate> pc;
    pp_max.pick(xrand, int_array, pc);

    TEST_EQUAL(pc.size(), 1)
    TEST_EQUAL(pc[0].pos, 10)
    TOLERANCE_RELATIVE(1.005);
    TEST_REAL_SIMILAR(pc[0].int_max, 0.07978845608028655)
    TOLERANCE_RELATIVE(1.02);
    TEST_REAL_SIMILAR(pc[0].mz_max, 10)
  }

  // Re-sample at every second point
  // Introduce a small int error
  {
    std::vector<PeakPickerMaxima::PeakCandidate> pc;
    pp_max.pick(mz_array, yrand, pc);

    TEST_EQUAL(pc.size(), 2)

    TEST_EQUAL(pc[0].pos, 8)
    TOLERANCE_RELATIVE(1.05);
    TEST_REAL_SIMILAR(pc[0].int_max, 0.07837041224348473) // yrand[8]
    TEST_EQUAL(fabs(pc[0].mz_max - 8) < 1, true)

    TEST_EQUAL(pc[1].pos, 11)
    TOLERANCE_RELATIVE(1.05);
    TEST_REAL_SIMILAR(pc[1].int_max, 0.07886053568902987) // yrand[10]
    TEST_EQUAL(fabs(pc[1].mz_max - 10) < 1, true)
  }

  TOLERANCE_RELATIVE(1.00001);
}
END_SECTION

START_SECTION([EXTRA](pick multiple peaks))
{
  PeakPickerMaxima pp_max(0.0);
  std::vector<double> mz_array, int_array;
  generateTestData(mz_array, int_array, 0, 1.0);
  for (Size i = 20; i < 25; i++)
  {
    mz_array.push_back(i);
    int_array.push_back( 0.020);
  }
  generateTestData(mz_array, int_array, 25, 1.0);
  for (Size i = 45; i < 50; i++)
  {
    mz_array.push_back(i);
    int_array.push_back( 0.020);
  }

  // Test multiple Gaussian function
  {
    std::vector<PeakPickerMaxima::PeakCandidate> pc;
    pp_max.pick(mz_array, int_array, pc);

    TEST_EQUAL(pc.size(), 2)
    TEST_EQUAL(pc[0].pos, 10)
    TEST_REAL_SIMILAR(pc[0].int_max, 0.07978845608028655)
    TEST_REAL_SIMILAR(pc[0].mz_max, 10)
    TEST_EQUAL(pc[0].left_boundary, 0)
    TEST_EQUAL(pc[0].right_boundary, 19)

    TEST_EQUAL(pc[1].pos, 35)
    TEST_REAL_SIMILAR(pc[1].int_max, 0.07978845608028655)
    TEST_REAL_SIMILAR(pc[1].mz_max, 35)
    TEST_EQUAL(pc[1].left_boundary, 25)
    TEST_EQUAL(pc[1].right_boundary, 44)
  }
}
END_SECTION

START_SECTION([EXTRA](pick multiple peaks SN))
{
  // Since S/N always returns a value > 1, we have to multiply our intensities
  // by a factor of 100.
  PeakPickerMaxima pp_max(1.0);
  std::vector<double> mz_array, int_array;
  generateTestData(mz_array, int_array, 0, 100.0);
  for (Size i = 20; i < 25; i++)
  {
    mz_array.push_back(i);
    int_array.push_back( 0.015*100);
  }
  generateTestData(mz_array, int_array, 25, 100.0);
  for (Size i = 45; i < 50; i++)
  {
    mz_array.push_back(i);
    int_array.push_back( 0.015*100);
  }

  // Test multiple Gaussian function with Signal to Noise
  {
    std::vector<PeakPickerMaxima::PeakCandidate> pc;
    pp_max.pick(mz_array, int_array, pc);

    TEST_EQUAL(pc.size(), 2)
    TEST_EQUAL(pc[0].pos, 10)
    TEST_REAL_SIMILAR(pc[0].int_max, 0.07978845608028655*100)
    TEST_REAL_SIMILAR(pc[0].mz_max, 10)
    TEST_EQUAL(pc[0].left_boundary, 1)
    TEST_EQUAL(pc[0].right_boundary, 17)

    TEST_EQUAL(pc[1].pos, 35)
    TEST_REAL_SIMILAR(pc[1].int_max, 0.07978845608028655*100)
    TEST_REAL_SIMILAR(pc[1].mz_max, 35)
    TEST_EQUAL(pc[1].left_boundary, 26)
    TEST_EQUAL(pc[1].right_boundary, 42)
  }
}
END_SECTION

TOLERANCE_RELATIVE(1.00001);
PeakMap input, output;

/////////////////////////
// ORBITRAP data tests //
/////////////////////////

// load Orbitrap input data
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_ppmax.mzML"),input);

////////////////////////////////////////////
// ORBITRAP test 1 (w/o noise estimation) //
////////////////////////////////////////////

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn0_out.mzML"),output);

START_SECTION([EXTRA](pick))
{
  PeakPickerMaxima pp_max(0.0);
  std::vector<PeakPickerMaxima::PeakCandidate> pc = ppmax_pick(input[0], pp_max);

  TOLERANCE_RELATIVE(1.1);
  // Check first scan
  TEST_EQUAL(pc.size(), 679)
  TEST_EQUAL(pc.size(), output[0].size())
  for (Size peak_idx = 0; peak_idx < output[0].size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(pc[peak_idx].mz_max, output[0][peak_idx].getMZ())
    TEST_REAL_SIMILAR(pc[peak_idx].int_max, output[0][peak_idx].getIntensity())
  }

  // Check all scans
  for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
  {
    pc = ppmax_pick(input[scan_idx], pp_max);
    TEST_EQUAL(output[scan_idx].size(), pc.size());
    for (Size peak_idx = 0; peak_idx < pc.size(); ++peak_idx)
    {
      TEST_REAL_SIMILAR(pc[peak_idx].mz_max, output[scan_idx][peak_idx].getMZ())
      TEST_REAL_SIMILAR(pc[peak_idx].int_max, output[scan_idx][peak_idx].getIntensity())
    }
  }
}
END_SECTION

/////////////////////////////////////////
// ORBITRAP test 2 (signal-to-noise 4) //
/////////////////////////////////////////

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn4_out_ppmax.mzML"),output);
START_SECTION([EXTRA](pick))
{
  PeakPickerMaxima pp_max(4.0);
  std::vector<PeakPickerMaxima::PeakCandidate> pc = ppmax_pick(input[0], pp_max);
  TOLERANCE_RELATIVE(1.05);

  // Check first scan
  TEST_EQUAL(output[0].size(), pc.size())
  for (Size peak_idx = 0; peak_idx < output[0].size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(pc[peak_idx].mz_max, output[0][peak_idx].getMZ())
    TEST_REAL_SIMILAR(pc[peak_idx].int_max, output[0][peak_idx].getIntensity())
  }

  // Check all scans
  for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
  {
    pc = ppmax_pick(input[scan_idx], pp_max);
    TEST_EQUAL(output[scan_idx].size(), pc.size())
    for (Size peak_idx = 0; peak_idx < output[scan_idx].size(); ++peak_idx)
    {
      TEST_REAL_SIMILAR(pc[peak_idx].mz_max, output[scan_idx][peak_idx].getMZ())
      TEST_REAL_SIMILAR(pc[peak_idx].int_max, output[scan_idx][peak_idx].getIntensity())
    }
  }
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/////////////////////////
// FTICR-MS data tests //
/////////////////////////

// load FTMS input data
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_ppmax.mzML"),input);

////////////////////////////////////////////
// FTICR-MS test 1 (w/o noise estimation) //
////////////////////////////////////////////

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn0_out.mzML"),output);
START_SECTION([EXTRA](pick))
{
  PeakPickerMaxima pp_max(0.0);
  std::vector<PeakPickerMaxima::PeakCandidate> pc = ppmax_pick(input[0], pp_max);

  // Check first scan
  TEST_EQUAL(pc.size(), 9359)
  TEST_EQUAL(output[0].size(), pc.size())
  int unequal_tests = 0;
  for (Size peak_idx = 0; peak_idx < output[0].size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(pc[peak_idx].mz_max, output[0][peak_idx].getMZ())
    if (fabs(pc[peak_idx].int_max - output[0][peak_idx].getIntensity())/output[0][peak_idx].getIntensity() > 0.05)
    {
      unequal_tests++;
    }
  }
  TEST_EQUAL(unequal_tests, 0)

  // Check all scans
  for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
  {
    pc = ppmax_pick(input[scan_idx], pp_max);
    TEST_EQUAL(output[scan_idx].size(), pc.size())
    for (Size peak_idx = 0; peak_idx < output[scan_idx].size(); ++peak_idx)
    {
      TEST_REAL_SIMILAR(pc[peak_idx].mz_max, output[scan_idx][peak_idx].getMZ())
      if (fabs(pc[peak_idx].int_max - output[scan_idx][peak_idx].getIntensity())/output[scan_idx][peak_idx].getIntensity() > 0.05)
      {
        unequal_tests++;
      }
    }
  }
  TEST_EQUAL(unequal_tests, 0)
}
END_SECTION

output.clear(true);

/////////////////////////////////////////
// FTICR-MS test 2 (signal-to-noise 4) //
/////////////////////////////////////////

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn4_out_ppmax.mzML"),output);

START_SECTION([EXTRA](template <typename PeakType> void pick(const MSSpectrum& input, MSSpectrum& output)))
{
  // With the new S/N the meaning of the noise value is slightly different:
  //  instead of the mean of the bin where the median can be found it is now
  //  the actual median. This new noise estimation seems to be generally lower
  //  than computed with the old method.
  //  -> to compensate we have to chose a _higher_ S/N cutoff 

  // Set the tolerance to 0.005 % 
  TOLERANCE_RELATIVE(1.000005);

  {
    PeakPickerMaxima pp_max(5.7);
    std::vector<PeakPickerMaxima::PeakCandidate> pc = ppmax_pick(input[0], pp_max);

    TEST_EQUAL(pc.size(), output[0].size() );
    int unequal_tests = 0;
    for (Size i = 0; i < pc.size(); i++)
    {
      TEST_REAL_SIMILAR(pc[i].mz_max,  output[0][i].getMZ())
      if (fabs(pc[i].int_max - output[0][i].getIntensity())/output[0][i].getIntensity() > 0.05)
      {
        unequal_tests++;
      }
    }
    TEST_EQUAL(unequal_tests, 0);
  }

  {
    PeakPickerMaxima pp_max(6.93);
    std::vector<PeakPickerMaxima::PeakCandidate> pc = ppmax_pick(input[1], pp_max);

    int unequal_tests = 0;
    TEST_EQUAL(pc.size(), output[1].size() );
    for (Size i = 0; i < pc.size(); i++)
    {
      TEST_REAL_SIMILAR(pc[i].mz_max,  output[1][i].getMZ())
      if (fabs(pc[i].int_max - output[1][i].getIntensity())/output[1][i].getIntensity() > 0.05)
      {
        unequal_tests++;
      }
    }
    TEST_EQUAL(unequal_tests, 0);
  }
}
END_SECTION

END_TEST

