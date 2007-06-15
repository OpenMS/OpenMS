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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestSignalToNoiseEstimator
  : public SignalToNoiseEstimator< >
{
  public:
  TestSignalToNoiseEstimator() 
    : SignalToNoiseEstimator< >()
  {
  }
  
  TestSignalToNoiseEstimator(const TestSignalToNoiseEstimator& bpf) 
  : SignalToNoiseEstimator< >(bpf)
  { 
  }
  
  TestSignalToNoiseEstimator& operator=(const TestSignalToNoiseEstimator& bpf)
  {
    if (&bpf==this) return *this;
    
    SignalToNoiseEstimator< >::operator=(bpf);
    
    return *this;
  }
  
  protected:
  
  virtual void computeSTN_(const PeakIterator& scan_first_, const PeakIterator& scan_last_) 
      throw(Exception::InvalidValue)
  {
    if (scan_first_ == scan_last_)
    {
      std::cout << "bla";
    }
    // do nothing here...
  }
  virtual DoubleReal getMaxIntensity() const { return 1;}
  virtual void setMaxIntensity(DoubleReal max_intensity)  { max_intensity = 1;}

  virtual DoubleReal getAutoMaxStdevFactor() const  { return 1;};
  virtual void setAutoMaxStdevFactor(DoubleReal value) { value = 1;}

  virtual DoubleReal getAutoMaxPercentile() const  { return 1;};
  virtual void setAutoMaxPercentile(DoubleReal value) { value = 1;}

  virtual inline Int getAutoMode() const  { return 1;};
  virtual void setAutoMode(Int value) { value = 1;}

  virtual DoubleReal getWinLen() const  { return 1;};
  virtual void setWinLen(DoubleReal win_len) { win_len = 1;}

  virtual Int getBinCount() const  { return 1;};
  virtual void setBinCount(Int bin_count) { bin_count = 1;}

  virtual Int getMinReqElements() const  { return 1;};
  virtual void setMinReqElements(Int min_required_elements) { min_required_elements = 1;}

  virtual DoubleReal getNoiseForEmtpyWindow() const  { return 1;};
  virtual void setNoiseForEmtpyWindow(DoubleReal noise_for_empty_window) { noise_for_empty_window = 1;} 
  

};

START_TEST(SignalToNoiseEstimator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSignalToNoiseEstimator* ptr = 0;
CHECK((SignalToNoiseEstimator()))
        ptr = new TestSignalToNoiseEstimator();
        TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((const PeakIterator& getFirstDataPoint() const))
  // done above
RESULT

CHECK((void setFirstDataPoint(const PeakIterator &first)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.setFirstDataPoint(spec.begin());
  
  TEST_EQUAL(sne.getFirstDataPoint() == spec.begin(), true)
RESULT

CHECK((const PeakIterator& getLastDataPoint() const))
  // done above
RESULT

CHECK((void setLastDataPoint(const PeakIterator &last)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.setLastDataPoint(spec.begin());
  
  TEST_EQUAL(sne.getLastDataPoint() == spec.begin(), true)
RESULT

CHECK((virtual DoubleReal getMaxIntensity() const =0))
{
  // abstract ...
}
RESULT

CHECK((virtual void setMaxIntensity(DoubleReal max_intensity)=0))
{
  // abstract ...
}
RESULT

CHECK((virtual DoubleReal getAutoMaxStdevFactor() const =0))
{
  // abstract ...
}
RESULT

CHECK((virtual void setAutoMaxStdevFactor(DoubleReal value)=0))
{
  // abstract ...
}
RESULT

CHECK((virtual DoubleReal getAutoMaxPercentile() const =0))
{
  // abstract ...
}
RESULT

CHECK((virtual void setAutoMaxPercentile(DoubleReal value)=0))
{
  // abstract ...
}
RESULT

CHECK((virtual Int getAutoMode() const =0))
{
  // abstract ...
}
RESULT

CHECK((virtual void setAutoMode(Int auto_mode)=0))
{
  // abstract ...
}
RESULT

CHECK((virtual DoubleReal getWinLen() const =0))
{
  // abstract ...
}
RESULT

CHECK((virtual void setWinLen(DoubleReal win_len)=0))
{
  // abstract ...
}
RESULT

CHECK((virtual Int getBinCount() const =0))
{
  // abstract ...
}
RESULT

CHECK((virtual void setBinCount(Int bin_count)=0))
{
  // abstract ...
}
RESULT

CHECK((virtual Int getMinReqElements() const =0))
{
  // abstract ...
}
RESULT

CHECK((virtual void setMinReqElements(Int min_required_elements)=0))
{
  // abstract ...
}
RESULT

CHECK((virtual DoubleReal getNoiseForEmtpyWindow() const =0))
{
  // abstract ...
}
RESULT

CHECK((virtual void setNoiseForEmtpyWindow(DoubleReal noise_for_empty_window)=0))
{
  // abstract ...
}
RESULT


CHECK((SignalToNoiseEstimator(const SignalToNoiseEstimator &source)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec.begin(), spec.end());
  
  
  TestSignalToNoiseEstimator sne_copy(sne);

  TEST_REAL_EQUAL(sne_copy.getFirstDataPoint() == spec.begin(), true)
  TEST_REAL_EQUAL(sne_copy.getLastDataPoint() == spec.end(), true)
RESULT

CHECK((SignalToNoiseEstimator& operator=(const SignalToNoiseEstimator &source)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec.begin(), spec.end());
  
  
  TestSignalToNoiseEstimator sne_copy;
  sne_copy = sne;
  
  TEST_REAL_EQUAL(sne_copy.getFirstDataPoint() == spec.begin(), true)
  TEST_REAL_EQUAL(sne_copy.getLastDataPoint() == spec.end(), true)
RESULT


CHECK((virtual ~SignalToNoiseEstimator()))
        delete ptr;
RESULT


CHECK((virtual void init(const PeakIterator& it_begin, const PeakIterator& it_end)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec.begin(), spec.end());
  
  TEST_EQUAL(sne.getFirstDataPoint() == spec.begin(), true)
  TEST_EQUAL(sne.getLastDataPoint() == spec.end(), true)
RESULT

CHECK((virtual void init(const Container& c)))
  TestSignalToNoiseEstimator sne;
  MSSpectrum<> spec;
  sne.init(spec);
  
  TEST_EQUAL(sne.getFirstDataPoint() == spec.begin(), true)
  TEST_EQUAL(sne.getLastDataPoint() == spec.end(), true)
RESULT

CHECK((virtual double getSignalToNoise(const PeakIterator& data_point)))
  // hard to do without implementing computeSTN_ properly
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


