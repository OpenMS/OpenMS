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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

//uncomment if the reference files should be re-written
//(only do this if you are sure that the PeakPickerHiRes is working correctly)
//#define WRITE_REF_FILES

START_TEST(PeakPickerHiRes, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerHiRes* ptr = nullptr;
PeakPickerHiRes* nullPointer = nullptr;
START_SECTION((PeakPickerHiRes()))
  ptr = new PeakPickerHiRes();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeakPickerHiRes()))
  delete ptr;
END_SECTION

PeakPickerHiRes pp_hires;
Param param;

PeakMap input, output;

/////////////////////////
// ORBITRAP data tests //
/////////////////////////


// load Orbitrap input data
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzML"),input);

/////////////////////////////////////////
// ORBITRAP test 1 (signal-to-noise 1) //
/////////////////////////////////////////

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn1_out.mzML"),output);

//set data type (this is not stored correctly in mzData)
for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
{
  output[scan_idx].setType(SpectrumSettings::CENTROID);
}

// PeakPickerHiRes config
param.setValue("signal_to_noise", 1.0);
pp_hires.setParameters(param);

START_SECTION((template <typename PeakType> void pick(const MSSpectrum& input, MSSpectrum& output) const))
  MSSpectrum tmp_spec;
  pp_hires.pick(input[0], tmp_spec);
#ifdef WRITE_REF_FILES
  PeakMap tmp_exp = input;
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
  {
    pp_hires.pick(input[scan_idx],tmp_spec);
    tmp_exp[scan_idx] = tmp_spec;
  }
  MzMLFile().store("./PeakPickerHiRes_orbitrap_sn1_out.mzML", tmp_exp);
#endif

  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
  }
END_SECTION

START_SECTION((template <typename PeakType> void pick(const MSSpectrum& input, MSSpectrum& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = true) const))
  MSSpectrum tmp_spec;
  std::vector<PeakPickerHiRes::PeakBoundary> tmp_boundaries;
  pp_hires.pick(input[0], tmp_spec, tmp_boundaries);
#ifdef WRITE_REF_FILES
  PeakMap tmp_exp = input;
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
  {
    pp_hires.pick(input[scan_idx],tmp_spec);
    tmp_exp[scan_idx] = tmp_spec;
  }
  MzMLFile().store("./PeakPickerHiRes_orbitrap_sn1_out.mzML", tmp_exp);
#endif

  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
  }

  TEST_REAL_SIMILAR(tmp_boundaries[25].mz_min, 359.728698730469)
  TEST_REAL_SIMILAR(tmp_boundaries[25].mz_max, 359.736419677734)
  TEST_REAL_SIMILAR(tmp_boundaries[26].mz_min, 360.155609130859)
  TEST_REAL_SIMILAR(tmp_boundaries[26].mz_max, 360.173675537109)

END_SECTION

START_SECTION([EXTRA](template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output)))
  // does the same as pick method for spectra
  NOT_TESTABLE
END_SECTION

START_SECTION([EXTRA](template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output, std::vector<std::vector<PeakBoundary> >& boundaries_spec, std::vector<std::vector<PeakBoundary> >& boundaries_chrom)))
  // does the same as pick method for spectra
  NOT_TESTABLE
END_SECTION

START_SECTION((template <typename PeakType, typename ChromatogramPeakT> void pickExperiment(const MSExperiment<PeakType, ChromatogramPeakT>& input, MSExperiment<PeakType, ChromatogramPeakT>& output) const))
  PeakMap tmp_exp;
  pp_hires.pickExperiment(input,tmp_exp);

  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
  {
    for (Size peak_idx = 0; peak_idx < tmp_exp[scan_idx].size(); ++peak_idx)
    {
      TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getMZ(), output[scan_idx][peak_idx].getMZ())
      TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getIntensity(), output[scan_idx][peak_idx].getIntensity())
    }
  }
END_SECTION

output.clear(true);

///////////////////////////////////////////
//// ORBITRAP test 2 (signal-to-noise 4) //
///////////////////////////////////////////


MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn4_out.mzML"),output);

//set data type (this is not stored correctly in mzData)
for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
{
  output[scan_idx].setType(SpectrumSettings::CENTROID);
}

//set up PeakPicker
param.setValue("signal_to_noise", 4.0);
pp_hires.setParameters(param);

START_SECTION([EXTRA](template <typename PeakType> void pick(const MSSpectrum& input, MSSpectrum& output)))
  MSSpectrum tmp_spec;
  pp_hires.pick(input[0],tmp_spec);
#ifdef WRITE_REF_FILES
  PeakMap tmp_exp = input;
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
  {
    pp_hires.pick(input[scan_idx],tmp_spec);
    tmp_exp[scan_idx] = tmp_spec;
  }
  MzMLFile().store("./PeakPickerHiRes_orbitrap_sn4_out.mzML", tmp_exp);
#endif

  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
  }
END_SECTION

START_SECTION([EXTRA](template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output)))
  PeakMap tmp_exp;
  pp_hires.pickExperiment(input,tmp_exp);

  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
  {
    for (Size peak_idx = 0; peak_idx < tmp_exp[scan_idx].size(); ++peak_idx)
    {
      TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getMZ(), output[scan_idx][peak_idx].getMZ())
      TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getIntensity(), output[scan_idx][peak_idx].getIntensity())
    }
  }
END_SECTION

output.clear(true);
input.clear(true);
//
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

/////////////////////////
// FTICR-MS data tests //
/////////////////////////


// load FTMS input data
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms.mzML"),input);

////////////////////////////////////////////////
//// FTICR-MS test 1 (signal-to-noise 1) //
////////////////////////////////////////////////

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn1_out.mzML"),output);

//set data type (this is not stored correctly in mzML)
for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
{
  output[scan_idx].setType(SpectrumSettings::CENTROID);
}

// PeakPickerHiRes config
param.setValue("signal_to_noise", 1.0);
pp_hires.setParameters(param);

START_SECTION([EXTRA](template <typename PeakType> void pick(const MSSpectrum& input, MSSpectrum& output)))
  MSSpectrum tmp_spec;
  pp_hires.pick(input[0],tmp_spec);
#ifdef WRITE_REF_FILES
  PeakMap tmp_exp = input;
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
  {
    pp_hires.pick(input[scan_idx],tmp_spec);
    tmp_exp[scan_idx] = tmp_spec;
  }
  MzMLFile().store("./PeakPickerHiRes_ftms_sn1_out.mzML", tmp_exp);
#endif
  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
  }
END_SECTION

output.clear(true);

/////////////////////////////////////////
// FTICR-MS test 2 (signal-to-noise 4) //
/////////////////////////////////////////


MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn4_out.mzML"),output);

//set data type (this is not stored correctly in mzML)
for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
{
  output[scan_idx].setType(SpectrumSettings::CENTROID);
}

//set up PeakPicker
param.setValue("signal_to_noise", 4.0);
pp_hires.setParameters(param);

START_SECTION([EXTRA](template <typename PeakType> void pick(const MSSpectrum& input, MSSpectrum& output)))
  MSSpectrum tmp_spec;
  pp_hires.pick(input[0],tmp_spec);
#ifdef WRITE_REF_FILES
  PeakMap tmp_exp = input;
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
  {
    pp_hires.pick(input[scan_idx],tmp_spec);
    tmp_exp[scan_idx] = tmp_spec;
  }
  MzMLFile().store("./PeakPickerHiRes_ftms_sn4_out.mzML", tmp_exp);
#endif

  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
  {
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
    TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
  }
END_SECTION

START_SECTION([EXTRA](template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output)))
  PeakMap tmp_exp;
  pp_hires.pickExperiment(input,tmp_exp);

  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
  {
    for (Size peak_idx = 0; peak_idx < tmp_exp[scan_idx].size(); ++peak_idx)
    {
      TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getMZ(), output[scan_idx][peak_idx].getMZ())
      TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getIntensity(), output[scan_idx][peak_idx].getIntensity())
    }
  }
END_SECTION

output.clear(true);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] test spectrum level selection)

  PeakMap inSpecSelection;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_spectrum_selection.mzML"), inSpecSelection);

  Param pp_hires_param;
  PeakPickerHiRes pp_spec_select;

  // pick only ms2
  PeakMap outMs2Only;
  pp_hires_param.setValue("ms_levels", ListUtils::create<Int>("2"));
  pp_spec_select.setParameters(pp_hires_param);

  pp_spec_select.pickExperiment(inSpecSelection, outMs2Only);

  ABORT_IF(inSpecSelection.size() != outMs2Only.size())
  for(Size i = 0; i < outMs2Only.size(); ++i)
  {
    if (outMs2Only[i].getMSLevel() == 2)
    {
      TEST_NOT_EQUAL(inSpecSelection[i], outMs2Only[i])
    }
    else
    {
      TEST_EQUAL(inSpecSelection[i], outMs2Only[i])
    }
  }

  // pick only ms1
  PeakMap outMs1Only;
  pp_hires_param.setValue("ms_levels", ListUtils::create<Int>("1"));
  pp_spec_select.setParameters(pp_hires_param);

  pp_spec_select.pickExperiment(inSpecSelection, outMs1Only);

  ABORT_IF(inSpecSelection.size() != outMs1Only.size())
  for(Size i = 0; i < outMs2Only.size(); ++i)
  {
    if (outMs2Only[i].getMSLevel() == 1)
    {
      TEST_NOT_EQUAL(inSpecSelection[i], outMs1Only[i])
    }
    else
    {
      TEST_EQUAL(inSpecSelection[i], outMs1Only[i])
    }
  }

  // pick ms1 and ms2
  PeakMap outMs1And2;
  pp_hires_param.setValue("ms_levels", ListUtils::create<Int>("1,2"));
  pp_spec_select.setParameters(pp_hires_param);

  pp_spec_select.pickExperiment(inSpecSelection, outMs1And2);

  ABORT_IF(inSpecSelection.size() != outMs2Only.size())
  for(Size i = 0; i < outMs2Only.size(); ++i)
  {
    if (outMs1And2[i].getMSLevel() == 2 || outMs1And2[i].getMSLevel() == 1)
    {
      TEST_NOT_EQUAL(inSpecSelection[i], outMs1And2[i])
    }
  }
END_SECTION

//////////////////////////////////////////////
// check peak boundaries on simulation data //
//////////////////////////////////////////////

// load input data
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_simulation.mzML"),input);

//set params
param.setValue("signal_to_noise", 0.0);
param.setValue("missing", 1);
param.setValue("spacing_difference_gap", 4.0);
pp_hires.setParameters(param);

START_SECTION(void pick(const MSSpectrum& input, MSSpectrum& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = true) const)
    PeakMap tmp_picked;
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > tmp_boundaries_s; // peak boundaries for spectra
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > tmp_boundaries_c; // peak boundaries for chromatograms

    pp_hires.pickExperiment(input, tmp_picked, tmp_boundaries_s, tmp_boundaries_c);

    TEST_EQUAL(tmp_picked[0].size(), 167);
    MSSpectrum::Iterator it_mz = tmp_picked.begin()->begin();
    vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundary = tmp_boundaries_s.begin()->begin();
    
    it_mz += 146;
    it_mz_boundary += 146;
    TEST_REAL_SIMILAR(it_mz->getMZ(),1141.57188829383);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_min,1141.51216791402);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_max,1141.63481354941);
    
    it_mz += 2;
    it_mz_boundary += 2;
    TEST_REAL_SIMILAR(it_mz->getMZ(),1142.57196823237);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_min,1142.50968574851);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_max,1142.6323313839);
    
    it_mz += 10;
    it_mz_boundary += 10;
    TEST_REAL_SIMILAR(it_mz->getMZ(),1178.08692219102);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_min,1178.02013862689);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_max,1178.14847787348);
    
    it_mz += 1;
    it_mz_boundary += 1;
    TEST_REAL_SIMILAR(it_mz->getMZ(),1178.58906411531);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_min,1178.5249396635);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_max,1178.6532789101);
    
END_SECTION

input.clear(true);
output.clear(true);

////////////////////////////////////////////
// check peak boundaries on orbitrap data //
////////////////////////////////////////////

// load input data
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzML"),input);

//set params
param.setValue("signal_to_noise", 0.0);
param.setValue("missing", 1);
param.setValue("spacing_difference_gap", 4.0);
pp_hires.setParameters(param);

START_SECTION(void pick(const MSSpectrum& input, MSSpectrum& output, std::vector<PeakBoundary>& boundaries, bool check_spacings = true) const)
    PeakMap tmp_picked;
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > tmp_boundaries_s; // peak boundaries for spectra
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > tmp_boundaries_c; // peak boundaries for chromatograms

    pp_hires.pickExperiment(input, tmp_picked, tmp_boundaries_s, tmp_boundaries_c);

    TEST_EQUAL(tmp_picked[0].size(), 82);
    MSSpectrum::Iterator it_mz = tmp_picked.begin()->begin();
    vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundary = tmp_boundaries_s.begin()->begin();
    
    it_mz += 14;
    it_mz_boundary += 14;
    TEST_REAL_SIMILAR(it_mz->getMZ(),355.070081088692);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_min,355.064544677734);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_max,355.078430175781);

    it_mz += 23;
    it_mz_boundary += 23;
    TEST_REAL_SIMILAR(it_mz->getMZ(),362.848715607077);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_min,362.844085693359);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_max,362.851928710938);

    it_mz += 17;
    it_mz_boundary += 17;
    TEST_REAL_SIMILAR(it_mz->getMZ(),370.210756298155);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_min,370.205871582031);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_max,370.215301513672); // Same as min of next peak.

    it_mz += 1;
    it_mz_boundary += 1;
    TEST_REAL_SIMILAR(it_mz->getMZ(),370.219596356153);
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_min,370.215301513672); // Same as max of previous peak.
    TEST_REAL_SIMILAR((*it_mz_boundary).mz_max,370.223358154297);
    
END_SECTION

END_TEST
