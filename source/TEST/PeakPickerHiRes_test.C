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
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerHiRes, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerHiRes* ptr = 0;
PeakPickerHiRes* nullPointer = 0;
START_SECTION((PeakPickerHiRes()))
  ptr = new PeakPickerHiRes();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PeakPickerHiRes()))
  delete ptr;
END_SECTION



PeakPickerHiRes pp_hires;
Param param;

MSExperiment<Peak1D> input, output;




/////////////////////////
// ORBITRAP data tests //
/////////////////////////


// load Orbitrap input data
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzML"),input);

////////////////////////////////////////////
// ORBITRAP test 1 (w/o noise estimation) //
////////////////////////////////////////////

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn0_out.mzML"),output);

//set data type (this is not stored correctly in mzData)
for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
{
	output[scan_idx].setType(SpectrumSettings::PEAKS);
}

// PeakPickerHiRes config
param.setValue("signal_to_noise",0.0);
pp_hires.setParameters(param);   

START_SECTION((template < typename PeakType > void pick(const MSSpectrum< PeakType > &input, MSSpectrum< PeakType > &output) const ))
  MSSpectrum<Peak1D> tmp_spec;
pp_hires.pick(input[0],tmp_spec);

// TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(output[0]), true)
  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
	{
		TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
		TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
	}
END_SECTION

START_SECTION((template < typename PeakType > void pickExperiment(const MSExperiment< PeakType > &input, MSExperiment< PeakType > &output) const ))
  MSExperiment<Peak1D> tmp_exp;
pp_hires.pickExperiment(input,tmp_exp);

TEST_EQUAL(tmp_exp.ExperimentalSettings::operator==(input), true)
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
	{
    // TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(output[scan_idx]), true)
		for (Size peak_idx = 0; peak_idx < tmp_exp[scan_idx].size(); ++peak_idx)
		{
			TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getMZ(), output[scan_idx][peak_idx].getMZ())
			TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getIntensity(), output[scan_idx][peak_idx].getIntensity())
		}
	}
END_SECTION

output.clear(true);



/////////////////////////////////////////
// ORBITRAP test 2 (signal-to-noise 4) //
/////////////////////////////////////////


MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn4_out.mzML"),output);

//set data type (this is not stored correctly in mzData)
for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
	{
		output[scan_idx].setType(SpectrumSettings::PEAKS);
	}

//set up PeakPicker
param.setValue("signal_to_noise",4.0);
pp_hires.setParameters(param);

START_SECTION([EXTRA](template <typename PeakType> void pick(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output)))
  MSSpectrum<Peak1D> tmp_spec;
pp_hires.pick(input[0],tmp_spec);

// TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(output[0]), true)
  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
		{
			TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
				TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
				}
END_SECTION

START_SECTION([EXTRA](template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output)))
  MSExperiment<Peak1D> tmp_exp;
pp_hires.pickExperiment(input,tmp_exp);

TEST_EQUAL(tmp_exp.ExperimentalSettings::operator==(input), true)
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
		{
    // TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(output[scan_idx]), true)
				for (Size peak_idx = 0; peak_idx < tmp_exp[scan_idx].size(); ++peak_idx)
					{
						TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getMZ(), output[scan_idx][peak_idx].getMZ())
							TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getIntensity(), output[scan_idx][peak_idx].getIntensity())
							}
		}
END_SECTION

output.clear(true);
input.clear(true);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

/////////////////////////
// FTICR-MS data tests //
/////////////////////////


// load FTMS input data
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms.mzML"),input);

////////////////////////////////////////////
// FTICR-MS test 1 (w/o noise estimation) //
////////////////////////////////////////////

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn0_out.mzML"),output);

//set data type (this is not stored correctly in mzML)
for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
	{
		output[scan_idx].setType(SpectrumSettings::PEAKS);
	}

// PeakPickerHiRes config
param.setValue("signal_to_noise",0.0);
pp_hires.setParameters(param);   

START_SECTION([EXTRA](template <typename PeakType> void pick(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output)))
  MSSpectrum<Peak1D> tmp_spec;
pp_hires.pick(input[0],tmp_spec);

// TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(output[0]), true)
  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
		{
			TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
				TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
				}
END_SECTION

START_SECTION([EXTRA](template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output)))
  MSExperiment<Peak1D> tmp_exp;
pp_hires.pickExperiment(input,tmp_exp);

TEST_EQUAL(tmp_exp.ExperimentalSettings::operator==(input), true)
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
		{
    // TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(output[scan_idx]), true)
				for (Size peak_idx = 0; peak_idx < tmp_exp[scan_idx].size(); ++peak_idx)
					{
						TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getMZ(), output[scan_idx][peak_idx].getMZ())
							TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getIntensity(), output[scan_idx][peak_idx].getIntensity())					
							}
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
		output[scan_idx].setType(SpectrumSettings::PEAKS);
	}

//set up PeakPicker
param.setValue("signal_to_noise",4.0);
pp_hires.setParameters(param);

START_SECTION([EXTRA](template <typename PeakType> void pick(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output)))
  MSSpectrum<Peak1D> tmp_spec;
pp_hires.pick(input[0],tmp_spec);

// TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(output[0]), true)
  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
		{
			TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
				TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
				}
END_SECTION

START_SECTION([EXTRA](template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output)))
  MSExperiment<Peak1D> tmp_exp;
pp_hires.pickExperiment(input,tmp_exp);

TEST_EQUAL(tmp_exp.ExperimentalSettings::operator==(input), true)
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
		{
    // TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(output[scan_idx]), true)
				for (Size peak_idx = 0; peak_idx < tmp_exp[scan_idx].size(); ++peak_idx)
					{
						TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getMZ(), output[scan_idx][peak_idx].getMZ())
							TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getIntensity(), output[scan_idx][peak_idx].getIntensity())
							}
		}
END_SECTION

output.clear(true);

/////////////////////////////////
// repeat test with RichPeak1D //
/////////////////////////////////

MSExperiment<RichPeak1D> inRich, outRich;

MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms.mzML"),inRich);
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn4_out.mzML"),outRich);

//set data type (this is not stored correctly in mzML)
for (Size scan_idx = 0; scan_idx < outRich.size(); ++scan_idx)
	{
		outRich[scan_idx].setType(SpectrumSettings::PEAKS);
	}

//set up PeakPicker
param.setValue("signal_to_noise",4.0);
pp_hires.setParameters(param);

START_SECTION([EXTRA](template <typename PeakType> void pick(const MSSpectrum<PeakType>& inRich, MSSpectrum<PeakType>& outRich)))
  MSSpectrum<RichPeak1D> tmp_spec;
pp_hires.pick(inRich[0],tmp_spec);

// TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(outRich[0]), true)
  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
		{
			TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), outRich[0][peak_idx].getMZ())
				TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), outRich[0][peak_idx].getIntensity())
				}
END_SECTION

START_SECTION([EXTRA](template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& inRich, MSExperiment<PeakType>& outRich)))
  MSExperiment<RichPeak1D> tmp_exp;
pp_hires.pickExperiment(inRich,tmp_exp);


TOLERANCE_RELATIVE(1e-4)
TEST_EQUAL(tmp_exp.ExperimentalSettings::operator==(inRich), true)
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
	{
    // TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(outRich[scan_idx]), true)
		for (Size peak_idx = 0; peak_idx < tmp_exp[scan_idx].size(); ++peak_idx)
		{
			TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getMZ(), outRich[scan_idx][peak_idx].getMZ())
			TEST_REAL_SIMILAR(tmp_exp[scan_idx][peak_idx].getIntensity(), outRich[scan_idx][peak_idx].getIntensity())
		}
	}
END_SECTION

inRich.clear(true);
outRich.clear(true);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST



