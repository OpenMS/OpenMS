// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakPickerHiRes, "$Id: PeakPickerHiRes_test.C 5049 2009-04-14 14:32:46Z ekenar $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakPickerHiRes* ptr = 0;
START_SECTION((PeakPickerHiRes()))
  ptr = new PeakPickerHiRes();
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~PeakPickerHiRes()))
  delete ptr;
END_SECTION

//load input and output data
MzDataFile mz_data_file;
MSExperiment<Peak1D> input, output;

// load input/output Orbitrap data
mz_data_file.load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap.mzData"),input);
mz_data_file.load(OPENMS_GET_TEST_DATA_PATH("PeakPickerHiRes_orbitrap_snt1_out.mzData"),output);
 

/////////////////////////////
// ORBITRAP test 1 (SNT=1) //
/////////////////////////////

//set data type (this is not stored correctly in mzData)
for (Size s=0; s<output.size(); ++s)
{
  output[s].setType(SpectrumSettings::PEAKS);
}

//set up PeakPicker
  PeakPickerHiRes pp;
  Param param;
  param.setValue("signal_to_noise",1.0);
  pp.setParameters(param);   

START_SECTION((void pick(const MSSpectrum<>& input, MSSpectrum<>& output)))
  MSSpectrum<> spec;
  pp.pick(input[0],spec);
  
  TEST_EQUAL(spec.SpectrumSettings::operator==(output[0]), true)
  for (Size p = 0; p < spec.size(); ++p)
  {
    TEST_REAL_SIMILAR(spec[p].getMZ(), output[0][p].getMZ())
    TEST_REAL_SIMILAR(spec[p].getIntensity(), output[0][p].getIntensity())
  }
END_SECTION

START_SECTION((void pickExperiment(const MSExperiment<>& input, MSExperiment<>& ouput)))
  MSExperiment<> exp;
  pp.pickExperiment(input,exp);

  TEST_EQUAL(exp.ExperimentalSettings::operator==(input), true)
  for (Size s=0; s<exp.size(); ++s)
  {
    TEST_EQUAL(exp[s].SpectrumSettings::operator==(output[s]), true)
    for (Size p=0; p<exp[s].size(); ++p)
    {
      TEST_REAL_SIMILAR(exp[s][p].getMZ(), output[s][p].getMZ())
      TEST_REAL_SIMILAR(exp[s][p].getIntensity(), output[s][p].getIntensity())
    }
  }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



