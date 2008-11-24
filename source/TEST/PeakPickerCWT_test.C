// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
mz_data_file.load("data/PeakPickerCWT.mzData",exp_raw);
START_SECTION((template<typename InputPeakContainer, typename OutputPeakContainer > void pick(const InputPeakContainer& input_peak_container, OutputPeakContainer& picked_peaks_container, int ms_level = 1)))
  MSSpectrum<> peaks;
  PeakPickerCWT pp;
    
  pp.pick(exp_raw[0],peaks);
  MSSpectrum<>::const_iterator it = peaks.begin();
END_SECTION

START_SECTION((template<typename InputPeakIterator, typename OutputPeakContainer  > void pick(InputPeakIterator first, InputPeakIterator last, OutputPeakContainer& picked_peak_container, int ms_level = 1)))
  MSSpectrum<> peaks;
  PeakPickerCWT pp;
  
  pp.pick(exp_raw[0].begin(),exp_raw[0].end(),peaks,1);
  MSSpectrum<>::const_iterator it = peaks.begin();
END_SECTION

Param param;
param.setValue("thresholds:peak_bound",1500.0);

START_SECTION((template<typename InputPeakType, typename OutputPeakType > void pickExperiment(const MSExperiment< InputPeakType >& ms_exp_raw, MSExperiment<OutputPeakType>& ms_exp_peaks)))
  MSExperiment<> peaks;
  PeakPickerCWT pp;
  pp.setParameters(param);
   
  pp.pickExperiment(exp_raw,peaks);
  TEST_EQUAL(peaks.size() == exp_raw.size(), true)
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 9)
  ExperimentalSettings e = peaks;
  TEST_EQUAL(e == exp_raw, true)
END_SECTION

		
START_SECTION((template<typename InputSpectrumIterator, typename OutputPeakType > void pickExperiment(InputSpectrumIterator first, InputSpectrumIterator last, MSExperiment<OutputPeakType>& ms_exp_peaks)))
	MSExperiment<Peak1D > exp_raw_ext;
  mz_data_file.load("data/PeakPickerCWT.mzData",exp_raw_ext);
	MSExperiment<> peaks;
  PeakPickerCWT pp;
  pp.setParameters(param);
   
  pp.pickExperiment(exp_raw_ext.begin(),exp_raw_ext.end(),peaks);
  TEST_EQUAL(peaks.size() == exp_raw_ext.size(), true)   
  TEST_EQUAL((peaks[0].size() + peaks[1].size()), 9)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



