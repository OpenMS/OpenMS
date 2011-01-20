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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzMLFile.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerCWT, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerCWT* ptr = 0;
START_SECTION((PeakPickerCWT()))
  ptr = new PeakPickerCWT();
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~PeakPickerCWT()))
  delete ptr;
END_SECTION

//load input and output data
MzMLFile mz_ml_file;
MSExperiment<Peak1D> input, output;
mz_ml_file.load(OPENMS_GET_TEST_DATA_PATH("PeakPickerCWT_test.mzML"),input);
mz_ml_file.load(OPENMS_GET_TEST_DATA_PATH("PeakPickerCWT_test_output.mzML"),output);
//set data type (this is not stored correctly in mzData)
for (Size s=0; s<output.size(); ++s)
{
  output[s].setType(SpectrumSettings::PEAKS);
}

//set up PeakPicker
  PeakPickerCWT pp;
  Param param;
  param.setValue("peak_width",0.15);
  param.setValue("signal_to_noise",3.);
  pp.setParameters(param);   

START_SECTION((void pick(const MSSpectrum<> &input, MSSpectrum<> &output)))
  MSSpectrum<> spec;
  pp.pick(input[0],spec);
  
// TEST_EQUAL(spec.SpectrumSettings::operator==(output[0]), true)-> are not equal as peak picking step is written to the spectrum settings
  for (Size p=0; p<spec.size(); ++p)
  {
    TEST_REAL_SIMILAR(spec[p].getMZ(), output[0][p].getMZ())
    TEST_REAL_SIMILAR(spec[p].getIntensity(), output[0][p].getIntensity())
  }
END_SECTION


START_SECTION((void pickExperiment(const MSExperiment<> &input, MSExperiment<> &output)))
  MSExperiment<> exp;
  pp.pickExperiment(input,exp);

  TEST_EQUAL(exp.ExperimentalSettings::operator==(input), true)
  for (Size s=0; s<exp.size(); ++s)
  {
		//    TEST_EQUAL(exp[s].SpectrumSettings::operator==(output[s]), true) -> are not equal as peak picking step is written to the spectrum settings
    for (Size p=0; p<exp[s].size(); ++p)
    {
      TEST_REAL_SIMILAR(exp[s][p].getMZ(), output[s][p].getMZ())
      TEST_REAL_SIMILAR(exp[s][p].getIntensity(), output[s][p].getIntensity())
    }
  }
END_SECTION

START_SECTION(DoubleReal estimatePeakWidth(const MSExperiment<>& input))
  DoubleReal peak_width = pp.estimatePeakWidth(input);
TEST_REAL_SIMILAR(peak_width,0.15)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



