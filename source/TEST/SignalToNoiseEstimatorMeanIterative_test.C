// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/DTAFile.h>

///////////////////////////
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMeanIterative.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

#define DEBUG_TEST
#undef DEBUG_TEST

START_TEST(SignalToNoiseEstimatorMeanIterative, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SignalToNoiseEstimatorMeanIterative< >* ptr = 0;
SignalToNoiseEstimatorMeanIterative< >* nullPointer = 0;
START_SECTION((SignalToNoiseEstimatorMeanIterative()))
        ptr = new SignalToNoiseEstimatorMeanIterative<>;
	TEST_NOT_EQUAL(ptr, nullPointer)
        SignalToNoiseEstimatorMeanIterative<> sne;
END_SECTION



START_SECTION((SignalToNoiseEstimatorMeanIterative& operator=(const SignalToNoiseEstimatorMeanIterative &source)))
  MSSpectrum < > raw_data;
  SignalToNoiseEstimatorMeanIterative<> sne;
	sne.init(raw_data);
  SignalToNoiseEstimatorMeanIterative<> sne2 = sne;
	NOT_TESTABLE
END_SECTION

START_SECTION((SignalToNoiseEstimatorMeanIterative(const SignalToNoiseEstimatorMeanIterative &source)))
  MSSpectrum < > raw_data;
  SignalToNoiseEstimatorMeanIterative<> sne;
	sne.init(raw_data);
  SignalToNoiseEstimatorMeanIterative<> sne2(sne);
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual ~SignalToNoiseEstimatorMeanIterative()))
  delete ptr;
END_SECTION


START_SECTION([EXTRA](virtual void init(const PeakIterator& it_begin, const PeakIterator& it_end)))

	TOLERANCE_ABSOLUTE(0.5)

  MSSpectrum < > raw_data;
  MSSpectrum< >::const_iterator it;
  DTAFile dta_file;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimator_test.dta"), raw_data);
  
    
  SignalToNoiseEstimatorMeanIterative<> sne;  
	Param p;
	p.setValue("win_len", 40.1);
	p.setValue("noise_for_empty_window", 2.0);
	p.setValue("min_required_elements", 10);
	sne.setParameters(p);
  sne.init(raw_data.begin(),raw_data.end());

  MSSpectrum < > stn_data;
  
#ifdef DEBUG_TEST
  MSSpectrum < > stn_data__;
#endif
  
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimatorMeanIterative_test.out"), stn_data);
  int i = 0;
  for (it=raw_data.begin();it!=raw_data.end(); ++it)
  {
    TEST_REAL_SIMILAR (stn_data[i].getIntensity(), sne.getSignalToNoise(it));
#ifdef DEBUG_TEST    
    Peak1D peak = (*it);
    peak.setIntensity(stn_data[i].getIntensity() / sne.getSignalToNoise(it));
    stn_data__.push_back(peak);
#endif    
    ++i;
  }
  
#ifdef DEBUG_TEST
  dta_file.store(OPENMS_GET_TEST_DATA_PATH("SignalToNoiseEstimatorMeanIterative_test.debug"), stn_data__);
#endif  
  
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


