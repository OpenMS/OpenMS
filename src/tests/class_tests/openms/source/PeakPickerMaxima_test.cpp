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
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerMaxima.h>
///////////////////////////

#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;


std::vector<PeakPickerMaxima::PeakCandidate> ppmax_pick(MSSpectrum<>& spec, PeakPickerMaxima& pp_max)
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

START_TEST(PeakPickerMaxima, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerMaxima* ptr = 0;
PeakPickerMaxima* nullPointer = 0;
START_SECTION((PeakPickerMaxima()))
  ptr = new PeakPickerMaxima(0,0,0);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeakPickerMaxima()))
  delete ptr;
END_SECTION

MSExperiment<Peak1D> input, output;

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

  // Check first scan
  TEST_EQUAL(pc.size(), 682)
  TEST_EQUAL(pc.size(), output[0].size())
  int unequal_tests = 0;
  for (Size peak_idx = 0; peak_idx < output[0].size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(pc[peak_idx].mz_max, output[0][peak_idx].getMZ())
    if (fabs(pc[peak_idx].int_max - output[0][peak_idx].getIntensity())/output[0][peak_idx].getIntensity() > 0.01)
    {
      unequal_tests++;
      // std::cout << pc[peak_idx].int_max << " != " <<  output[0][peak_idx].getIntensity() << std::endl;
    }
  }
  TEST_EQUAL(unequal_tests, 26)

  TOLERANCE_RELATIVE(1.05);
  // Check all scans
  for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
  {
    pc = ppmax_pick(input[scan_idx], pp_max);
    TEST_EQUAL(output[scan_idx].size(), pc.size())
    for (Size peak_idx = 0; peak_idx < output[scan_idx].size(); ++peak_idx)
    {
      TEST_REAL_SIMILAR(pc[peak_idx].mz_max, output[scan_idx][peak_idx].getMZ())

      if (fabs(pc[peak_idx].int_max - output[scan_idx][peak_idx].getIntensity())/output[scan_idx][peak_idx].getIntensity() > 0.01)
      {
        unequal_tests++;
      }
    }
  }
  TEST_EQUAL(unequal_tests, 97)
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
  TEST_EQUAL(unequal_tests, 54)

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
  TEST_EQUAL(unequal_tests, 148)
}
END_SECTION

output.clear(true);

/////////////////////////////////////////
// FTICR-MS test 2 (signal-to-noise 4) //
/////////////////////////////////////////

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn4_out_ppmax.mzML"),output);

START_SECTION([EXTRA](template <typename PeakType> void pick(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output)))
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
      if (fabs(pc[i].int_max - output[0][i].getIntensity())/output[0][i].getIntensity() > 0.01)
      {
        unequal_tests++;
      }
    }
    TEST_EQUAL(unequal_tests, 13);
  }

  {
    PeakPickerMaxima pp_max(6.93);
    std::vector<PeakPickerMaxima::PeakCandidate> pc = ppmax_pick(input[1], pp_max);

    int unequal_tests = 0;
    TEST_EQUAL(pc.size(), output[1].size() );
    for (Size i = 0; i < pc.size(); i++)
    {
      TEST_REAL_SIMILAR(pc[i].mz_max,  output[1][i].getMZ())
      if (fabs(pc[i].int_max - output[1][i].getIntensity())/output[1][i].getIntensity() > 0.01)
      {
        unequal_tests++;
      }
    }
    TEST_EQUAL(unequal_tests, 7);
  }
}
END_SECTION

END_TEST

