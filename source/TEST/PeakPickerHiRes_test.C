// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzDataFile.h>

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
MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzData"),input);

////////////////////////////////////////////
// ORBITRAP test 1 (w/o noise estimation) //
////////////////////////////////////////////

MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn0_out.mzData"),output);

//set data type (this is not stored correctly in mzData)
for (Size scan_idx = 0; scan_idx < output.size(); ++scan_idx)
	{
		output[scan_idx].setType(SpectrumSettings::PEAKS);
	}

// PeakPickerHiRes config
param.setValue("signal_to_noise",0.0);
pp_hires.setParameters(param);   

START_SECTION((template <typename PeakType> void pick(const MSSpectrum<PeakType>& input, MSSpectrum<PeakType>& output)))
  MSSpectrum<Peak1D> tmp_spec;
pp_hires.pick(input[0],tmp_spec);

TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(output[0]), true)
  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
		{
			TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), output[0][peak_idx].getMZ())
				TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), output[0][peak_idx].getIntensity())
				}
END_SECTION

START_SECTION((template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& input, MSExperiment<PeakType>& output)))
  MSExperiment<Peak1D> tmp_exp;
pp_hires.pickExperiment(input,tmp_exp);

TEST_EQUAL(tmp_exp.ExperimentalSettings::operator==(input), true)
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
		{
			TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(output[scan_idx]), true)
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


MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_sn4_out.mzData"),output);

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

TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(output[0]), true)
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
			TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(output[scan_idx]), true)
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
MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms.mzData"),input);

////////////////////////////////////////////
// FTICR-MS test 1 (w/o noise estimation) //
////////////////////////////////////////////

MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn0_out.mzData"),output);

//set data type (this is not stored correctly in mzData)
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

TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(output[0]), true)
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
			TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(output[scan_idx]), true)
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


MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn4_out.mzData"),output);

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

TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(output[0]), true)
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
			TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(output[scan_idx]), true)
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

MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms.mzData"),inRich);
MzDataFile().load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_ftms_sn4_out.mzData"),outRich);

//set data type (this is not stored correctly in mzData)
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

TEST_EQUAL(tmp_spec.SpectrumSettings::operator==(outRich[0]), true)
  for (Size peak_idx = 0; peak_idx < tmp_spec.size(); ++peak_idx)
		{
			TEST_REAL_SIMILAR(tmp_spec[peak_idx].getMZ(), outRich[0][peak_idx].getMZ())
				TEST_REAL_SIMILAR(tmp_spec[peak_idx].getIntensity(), outRich[0][peak_idx].getIntensity())
				}
END_SECTION

START_SECTION([EXTRA](template <typename PeakType> void pickExperiment(const MSExperiment<PeakType>& inRich, MSExperiment<PeakType>& outRich)))
  MSExperiment<RichPeak1D> tmp_exp;
pp_hires.pickExperiment(inRich,tmp_exp);

TEST_EQUAL(tmp_exp.ExperimentalSettings::operator==(inRich), true)
  for (Size scan_idx = 0; scan_idx < tmp_exp.size(); ++scan_idx)
		{
			TEST_EQUAL(tmp_exp[scan_idx].SpectrumSettings::operator==(outRich[scan_idx]), true)
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



