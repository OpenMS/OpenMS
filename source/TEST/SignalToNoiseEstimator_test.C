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
// $Maintainer: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

///////////////////////////
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestSignalToNoiseEstimator
  : public SignalToNoiseEstimator
{
  public:
  TestSignalToNoiseEstimator() 
    : SignalToNoiseEstimator()
  {
  }
  
  TestSignalToNoiseEstimator(const TestSignalToNoiseEstimator& bpf) 
  : SignalToNoiseEstimator(bpf)
  { 
  }
  
  TestSignalToNoiseEstimator& operator=(const TestSignalToNoiseEstimator& bpf)
  {
    if (&bpf==this) return *this;
    
    SignalToNoiseEstimator::operator=(bpf);
    
    return *this;
  }
  
  virtual double getSignalToNoise(PeakIterator data_point)
  {
    return data_point->getMZ();
  }

};

START_TEST(SignalToNoiseEstimator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSignalToNoiseEstimator* ptr = 0;
CHECK(SignalToNoiseEstimator())
        ptr = new TestSignalToNoiseEstimator();
        TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~SignalToNoiseEstimator())
        delete ptr;
RESULT

CHECK(virtual void init(PeakIterator it_begin, PeakIterator it_end))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec.begin(), spec.end());
  
  TEST_EQUAL(sne.getFirstDataPoint() == spec.begin(), true)
  TEST_EQUAL(sne.getLastDataPoint() == spec.end(), true)
RESULT

CHECK(const PeakIterator& getFirstDataPoint() const )
  // done above
RESULT

CHECK(void setFirstDataPoint(const PeakIterator &first))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.setFirstDataPoint(spec.begin());
  
  TEST_EQUAL(sne.getFirstDataPoint() == spec.begin(), true)
RESULT

CHECK(const PeakIterator& getLastDataPoint() const )
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.setLastDataPoint(spec.begin());
  
  TEST_EQUAL(sne.getLastDataPoint() == spec.begin(), true)
RESULT

CHECK(void setLastDataPoint(const PeakIterator &last))
  // done above
RESULT

CHECK(virtual double getSignalToNoise(PeakIterator data_point)=0)
  // just a base class... no implementation here
RESULT

CHECK(SignalToNoiseEstimator(const SignalToNoiseEstimator &source))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec.begin(), spec.end());
  
  
  TestSignalToNoiseEstimator sne_copy(sne);

  TEST_REAL_EQUAL(sne_copy.getFirstDataPoint() == spec.begin(), true)
  TEST_REAL_EQUAL(sne_copy.getLastDataPoint() == spec.end(), true)
RESULT

CHECK(SignalToNoiseEstimator& operator=(const SignalToNoiseEstimator &source))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec.begin(), spec.end());
  
  
  TestSignalToNoiseEstimator sne_copy;
  sne_copy = sne;
  
  TEST_REAL_EQUAL(sne_copy.getFirstDataPoint() == spec.begin(), true)
  TEST_REAL_EQUAL(sne_copy.getLastDataPoint() == spec.end(), true)
RESULT

CHECK(virtual ~SignalToNoiseEstimator())
  // ...
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


