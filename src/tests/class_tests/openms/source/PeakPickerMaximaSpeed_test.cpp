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

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

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

START_TEST(PeakPickerMaximaSpeed, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSExperiment<Peak1D> input;
//MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms.mzML"),input);
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_ppmax.mzML"),input);

START_SECTION(Speed_no_SN)
{
  // Set the tolerance to 0.005 % 
  TOLERANCE_RELATIVE(1.00005);

  MSSpectrum<Peak1D> tmp_spec;
  PeakPickerMaxima pp_max(0.0);
  PeakPickerHiRes pp_hires;
  Param param;
  param.setValue("signal_to_noise", 0.0);
  pp_hires.setParameters(param);   

  pp_hires.pick(input[0],tmp_spec);

  {
  int nr_peaks_picked = 0;
  clock_t begin = clock();
  for (Size i = 0; i < 1; i++)
  {
    pp_hires.pick(input[0],tmp_spec); nr_peaks_picked += tmp_spec.size();
    pp_hires.pick(input[1],tmp_spec); nr_peaks_picked += tmp_spec.size();
  }
  clock_t end = clock();
  std::cout << " Old Peakpicker time: " << double(end - begin) / CLOCKS_PER_SEC << " for " << nr_peaks_picked << " peaks." << std::endl;
  }

  int nr_peaks_picked = 0;
  clock_t begin = clock();
  // For this many peaks, it seems that the spline fitting also starts to
  // become a bottle neck (takes nearly 50% of the time)
  for (Size i = 0; i < 100; i++)
  {
    nr_peaks_picked += ppmax_pick(input[0], pp_max).size();
    nr_peaks_picked += ppmax_pick(input[1], pp_max).size();
  }
  clock_t end = clock();
  std::cout << " New Peakpicker time: " << double(end - begin) / CLOCKS_PER_SEC << " for " << nr_peaks_picked << " peaks." << std::endl;
}
END_SECTION

START_SECTION(Speed_SN)
{
  // Set the tolerance to 0.005 % 
  TOLERANCE_RELATIVE(1.00005);

  MSSpectrum<Peak1D> tmp_spec;
  PeakPickerMaxima pp_max(5.7);
  PeakPickerHiRes pp_hires;
  Param param;
  param.setValue("signal_to_noise", 4.0);
  pp_hires.setParameters(param);   

  pp_hires.pick(input[0],tmp_spec);

  {
  int nr_peaks_picked = 0;
  clock_t begin = clock();
  for (Size i = 0; i < 10; i++)
  {
    pp_hires.pick(input[0],tmp_spec); nr_peaks_picked += tmp_spec.size();
    pp_hires.pick(input[1],tmp_spec); nr_peaks_picked += tmp_spec.size();
  }
  clock_t end = clock();
  std::cout << " Old Peakpicker time: " << double(end - begin) / CLOCKS_PER_SEC << " for " << nr_peaks_picked << " peaks." << std::endl;
  }

  int nr_peaks_picked = 0;
  clock_t begin = clock();
  for (Size i = 0; i < 100; i++)
  {
    nr_peaks_picked += ppmax_pick(input[0], pp_max).size();
    nr_peaks_picked += ppmax_pick(input[1], pp_max).size();
  }
  clock_t end = clock();
  std::cout << " New Peakpicker time: " << double(end - begin) / CLOCKS_PER_SEC << " for " << nr_peaks_picked << " peaks." << std::endl;
}
END_SECTION

END_TEST

