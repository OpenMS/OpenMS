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
// $Maintainer: Eva Lange $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzDataFile.h>

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

MzDataFile mz_data_file;
MSExperiment<Peak1D > exp_raw;
mz_data_file.load(OPENMS_GET_TEST_DATA_PATH("PeakPickerCWT.mzData"),exp_raw);
START_SECTION((void pick(const MSSpectrum<>& input, MSSpectrum<>& ouput)))
  MSSpectrum<> peaks;
  PeakPickerCWT pp;
    
  pp.pick(exp_raw[0],peaks);
END_SECTION

Param param;
param.setValue("thresholds:peak_bound",1500.0);

START_SECTION((void pickExperiment(const MSExperiment<>& input, MSExperiment<>& ouput)))
  MSExperiment<> peaks;
  PeakPickerCWT pp;
  pp.setParameters(param);
   
  pp.pickExperiment(exp_raw,peaks);
  TEST_EQUAL(peaks.size() == exp_raw.size(), true)
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 9)
  ExperimentalSettings e = peaks;
  TEST_EQUAL(e == exp_raw, true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



