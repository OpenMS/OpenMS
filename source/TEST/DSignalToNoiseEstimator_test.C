// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/FILTERING/NOISEESTIMATION/DSignalToNoiseEstimator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef RawSpectrum::const_iterator PeakIteratorType;

class TestSignalToNoiseEstimator : public DSignalToNoiseEstimator<1, PeakIteratorType>
{
  public:
  TestSignalToNoiseEstimator() : DSignalToNoiseEstimator<1, RawSpectrum::const_iterator>(){}
  TestSignalToNoiseEstimator(const Param& param) : DSignalToNoiseEstimator<1, RawSpectrum::const_iterator>(param){}
  TestSignalToNoiseEstimator(const TestSignalToNoiseEstimator& bpf) : DSignalToNoiseEstimator<1, RawSpectrum::const_iterator>(bpf){}
  TestSignalToNoiseEstimator& operator=(const TestSignalToNoiseEstimator& bpf)
  {
     DSignalToNoiseEstimator<1, PeakIteratorType>::operator=(bpf);
     return *this;
  }
  virtual double getSignalToNoise(PeakIteratorType /* data_point*/)
  {
    return 0.;
  }

};

START_TEST(DSignalToNoiseEstimator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSignalToNoiseEstimator* ptr = 0;
CHECK((DSignalToNoiseEstimator()))
	ptr = new TestSignalToNoiseEstimator();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~DSignalToNoiseEstimator()))
	delete ptr;
RESULT

CHECK((DSignalToNoiseEstimator& operator=(const DSignalToNoiseEstimator& ne)))
  TestSignalToNoiseEstimator sne;
  sne.getMZdim()=-1;
  sne.getRTdim()=0;
  
  TestSignalToNoiseEstimator sne_copy;
  sne_copy = sne;
  
  TEST_REAL_EQUAL(sne_copy.getMZdim(),-1)
  TEST_REAL_EQUAL(sne_copy.getRTdim(),0)
RESULT

CHECK((DSignalToNoiseEstimator(const Param& parameters)))
  Param param;
  param.setValue("bla",3);
  TestSignalToNoiseEstimator sne(param);
  
  TEST_EQUAL(sne.getParam() == param, true)
RESULT

CHECK((DSignalToNoiseEstimator(const DSignalToNoiseEstimator& ne)))
  TestSignalToNoiseEstimator sne;
  sne.getMZdim()=-1;
  sne.getRTdim()=0;
  
  TestSignalToNoiseEstimator sne_copy(sne);
  
  TEST_REAL_EQUAL(sne_copy.getMZdim(),-1)
  TEST_REAL_EQUAL(sne_copy.getRTdim(),0)
RESULT

CHECK((PeakIterator& getFirstDataPoint()))
  TestSignalToNoiseEstimator sne;
  RawSpectrum spec;
  sne.getFirstDataPoint() = spec.begin();
  
  TEST_EQUAL(sne.getFirstDataPoint() == spec.begin(), true)
RESULT

CHECK((PeakIterator& getLastDataPoint()))
  TestSignalToNoiseEstimator sne;
  RawSpectrum spec;
  sne.getLastDataPoint() = spec.end();
  
  TEST_EQUAL(sne.getLastDataPoint() == spec.end(), true)
RESULT

CHECK((const Param& getParam() const))
  Param param;
  param.setValue("bla",3);
  TestSignalToNoiseEstimator sne(param);
  
  TEST_EQUAL(sne.getParam() == param, true)
RESULT

CHECK((const PeakIterator& getFirstDataPoint() const))
  TestSignalToNoiseEstimator sne;
  PeakIteratorType pit(0);
  
  TEST_EQUAL(sne.getFirstDataPoint() == pit, true)
RESULT

CHECK((const PeakIterator& getLastDataPoint() const))
  TestSignalToNoiseEstimator sne;
  PeakIteratorType pit(0);
  
  TEST_EQUAL(sne.getLastDataPoint() == pit, true)
RESULT

CHECK((const int getRTdim() const))
  TestSignalToNoiseEstimator sne;
  
  TEST_REAL_EQUAL(sne.getRTdim(),-1)
RESULT

CHECK((const int& getMZdim() const))
  TestSignalToNoiseEstimator sne;
  
  TEST_REAL_EQUAL(sne.getMZdim(),0)
RESULT

CHECK((double getSignalToNoise(PeakIterator data_point)))
 
RESULT

CHECK((int& getMZdim()))
  TestSignalToNoiseEstimator sne;
  sne.getMZdim()=-1;
  
  TEST_REAL_EQUAL(sne.getMZdim(),-1)
RESULT

CHECK((int& getRTdim()))
  TestSignalToNoiseEstimator sne;
  sne.getRTdim()=0;
  
  TEST_REAL_EQUAL(sne.getRTdim(),0)
RESULT

CHECK((void init(PeakIterator it_begin, PeakIterator it_end)))
  TestSignalToNoiseEstimator sne;
  RawSpectrum spec;
  sne.init(spec.begin(),spec.end());
  
  TEST_EQUAL(sne.getFirstDataPoint() == spec.begin(), true)
RESULT

CHECK((void setFirstDataPoint(const PeakIterator& first)))
  TestSignalToNoiseEstimator sne;
  RawSpectrum spec;
  sne.setFirstDataPoint(spec.begin());
  
  TEST_EQUAL(sne.getFirstDataPoint() == spec.begin(), true)
RESULT

CHECK((void setLastDataPoint(const PeakIterator& last)))
  TestSignalToNoiseEstimator sne;
  RawSpectrum spec;
  sne.setLastDataPoint(spec.end());
  
  TEST_EQUAL(sne.getLastDataPoint() == spec.end(), true)
RESULT

CHECK((void setMZdim(const int& mz_dim)))
  TestSignalToNoiseEstimator sne;
  sne.setMZdim(-1);
  
  TEST_REAL_EQUAL(sne.getMZdim(),-1)
RESULT

CHECK((void setParam(const Param& param)))
  Param param;
  param.setValue("bla",3);
  TestSignalToNoiseEstimator sne;
  sne.setParam(param);
  
  TEST_EQUAL(sne.getParam() == param, true)
RESULT

CHECK((void setRTdim(const int& rt_dim)))
  TestSignalToNoiseEstimator sne;
  sne.setRTdim(0);
  
  TEST_REAL_EQUAL(sne.getRTdim(),0)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



