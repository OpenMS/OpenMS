// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ContinuousWaveletTransform, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
class Test: public ContinuousWaveletTransform<1>
{
  public:
  Test() : ContinuousWaveletTransform<1>()
  {
  };

  Test(const Test& source): ContinuousWaveletTransform<1>(source)
  {
  };
		
  virtual ~Test()
  {		
  };
	
  Test& operator = (const Test& source)
  {
    if (&source != this)
    {
      ContinuousWaveletTransform<1>::operator=(source);
    }
    return *this;
  };
	
   virtual void transform(RawDataPointConstIterator /* begin_input */, RawDataPointConstIterator /* end_input */, float /*resolution*/)
   {
   };
};

Test* ptr = 0;
CHECK(Test())
	ptr = new Test();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~Test())
	delete ptr;
RESULT

CHECK(Test& operator=(const Test& cwt))
  Test cwt;
  
  cwt.getScale() = 0.2;
  cwt.getSpacing() = 0.01;
  cwt.getLeftPaddingIndex() = 3;
  cwt.getRightPaddingIndex() = 1;
  cwt.getSignalLength() = 2;
  cwt.getMzDim() = 1;
  
  Test cwt_copy;
  cwt_copy = cwt;
  
  TEST_REAL_EQUAL(cwt_copy.getScale(), 0.2)
  TEST_REAL_EQUAL(cwt_copy.getSpacing(), 0.01)
  TEST_REAL_EQUAL(cwt_copy.getLeftPaddingIndex(),3)
  TEST_REAL_EQUAL(cwt_copy.getRightPaddingIndex(),1)
  TEST_REAL_EQUAL(cwt_copy.getSignalLength(), 2)
  TEST_REAL_EQUAL(cwt_copy.getMzDim(),1)
RESULT

CHECK(Test(const Test& cwt))
  Test cwt;
  
  cwt.getScale() = 0.2;
  cwt.getSpacing() = 0.01;
  cwt.getLeftPaddingIndex() = 3;
  cwt.getRightPaddingIndex() = 1;
  cwt.getSignalLength() = 2;
  cwt.getMzDim() = 1;
  
  Test cwt_copy(cwt);
  
  TEST_REAL_EQUAL(cwt_copy.getScale(), 0.2)
  TEST_REAL_EQUAL(cwt_copy.getSpacing(), 0.01)
  TEST_REAL_EQUAL(cwt_copy.getLeftPaddingIndex(),3)
  TEST_REAL_EQUAL(cwt_copy.getRightPaddingIndex(),1)
  TEST_REAL_EQUAL(cwt_copy.getSignalLength(), 2)
  TEST_REAL_EQUAL(cwt_copy.getMzDim(),1)
RESULT

CHECK((DPeakArrayNonPolymorphic<1, DRawDataPoint<1>& getSignal()))
  Test cwt;
  
  DPeakArrayNonPolymorphic<1, DRawDataPoint<1> > signal(2);
  cwt.getSignal() = signal;
  
  TEST_EQUAL(cwt.getSignal() == signal, true)
RESULT

CHECK((const DPeakArrayNonPolymorphic<1, DRawDataPoint<1>& getSignal() const))
 Test cwt;
 
 TEST_EQUAL(cwt.getSignal().size() == 0, true)
RESULT

CHECK(const double& getScale() const)
  Test cwt;
  
  TEST_REAL_EQUAL(cwt.getScale(), 0)
RESULT

CHECK(const double& getSpacing() const)
  Test cwt;
  
  TEST_REAL_EQUAL(cwt.getSpacing(), 0)
RESULT

CHECK(const double& operator [] (const unsigned int i) const)
  Test cwt;
  
  DPeakArrayNonPolymorphic<1, DRawDataPoint<1> > signal(1);
  cwt.getSignal() = signal;
  
  Test const cwt_const(cwt);
  
  TEST_REAL_EQUAL(cwt_const[0],0)
RESULT

CHECK(const int& getLeftPaddingIndex() const)
  Test cwt;
  
  TEST_REAL_EQUAL(cwt.getLeftPaddingIndex(), 0)
RESULT

CHECK(const int& getRightPaddingIndex() const)
  Test cwt;
  
  TEST_REAL_EQUAL(cwt.getRightPaddingIndex(), 0)
RESULT

CHECK(const int& getSignalLength() const)
  Test cwt;
  
  TEST_REAL_EQUAL(cwt.getSignalLength(), 0)
RESULT

CHECK(const int& getSize() const)
  Test cwt;
  
  TEST_REAL_EQUAL(cwt.getSize(), 0)
RESULT

CHECK(const std::vector<double>& getWavelet() const)
  Test cwt;
  
  TEST_REAL_EQUAL(cwt.getScale(), 0)
RESULT

CHECK(const unsigned int& getMzDim() const)
  Test cwt;
  
  TEST_REAL_EQUAL(cwt.getScale(), 0)
RESULT

CHECK(double& getScale())
  Test cwt;
  cwt.getScale() = 0.2;
  
  TEST_REAL_EQUAL(cwt.getScale(), 0.2)
RESULT

CHECK(double& getSpacing())
  Test cwt;
  cwt.getSpacing() = 0.2;
  
  TEST_REAL_EQUAL(cwt.getSpacing(), 0.2)
RESULT

CHECK(double& operator [] (const unsigned int i))
  DPeakArrayNonPolymorphic<1, DRawDataPoint<1> > signal;
  DRawDataPoint<1> rp;
  rp.setIntensity(100);
  signal.push_back(rp);
  
  Test cwt;
  cwt.getSignal() = signal;
  
  TEST_EQUAL(cwt[0],100)
RESULT

CHECK(int& getLeftPaddingIndex())
  Test cwt;
  cwt.getLeftPaddingIndex() = 2;
  
  TEST_REAL_EQUAL(cwt.getLeftPaddingIndex(), 2)
RESULT

CHECK(int& getRightPaddingIndex())
  Test cwt;
  cwt.getRightPaddingIndex() = 2;
  
  TEST_REAL_EQUAL(cwt.getRightPaddingIndex(), 2)
RESULT

CHECK(int& getSignalLength())
  Test cwt;
  cwt.getSignalLength() = 2;
  
  TEST_REAL_EQUAL(cwt.getSignalLength(), 2)
RESULT

CHECK(std::vector<double>& getWavelet())
  vector<double> w(1);
  w[1] = 0.5;
  
  Test cwt;
  cwt.getWavelet() = w;
  
  TEST_EQUAL(cwt.getWavelet() == w, true)
RESULT

CHECK(unsigned int& getMzDim())
  Test cwt;
  cwt.getMzDim() = 1;
  
  TEST_REAL_EQUAL(cwt.getMzDim(), 1)
RESULT

CHECK((void init(double scale, double spacing, unsigned int MZ)))
  Test cwt;
  double scale = 0.2;
  double spacing = 2.3;
  unsigned int mz = 2;
  cwt.init(scale,spacing,mz);
  
  TEST_EQUAL(cwt.getMzDim(), mz)
  TEST_REAL_EQUAL(cwt.getSpacing(),spacing)
  TEST_REAL_EQUAL(cwt.getScale(),scale)
RESULT

CHECK(void setLeftPaddingIndex(const int end_left_padding))
  Test cwt;
  cwt.setLeftPaddingIndex(2);
  
  TEST_REAL_EQUAL(cwt.getLeftPaddingIndex(), 2)
RESULT

CHECK(void setMzDim(const unsigned int& mz_dim))
  Test cwt;
  cwt.setMzDim(2);
  
  TEST_EQUAL(cwt.getMzDim(), 2)
RESULT

CHECK(void setRightPaddingIndex(const int begin_right_padding))
  Test cwt;
  cwt.setRightPaddingIndex(2);
  
  TEST_EQUAL(cwt.getRightPaddingIndex(), 2)
RESULT

CHECK(void setScale(const double& scale))
  Test cwt;
  cwt.setScale(0.2);
  
  TEST_REAL_EQUAL(cwt.getScale(), 0.2)
RESULT

CHECK((void setSignal(const DPeakArrayNonPolymorphic<1, DRawDataPoint<1> >& signal)))
  Test cwt;
  
  DPeakArrayNonPolymorphic<1, DRawDataPoint<1> > signal(2);
  cwt.setSignal(signal);
  
  TEST_EQUAL(cwt.getSignal() == signal, true)
RESULT

CHECK(void setSignalLength(const int signal_length))
  Test cwt;
  cwt.setSignalLength(2);
  
  TEST_EQUAL(cwt.getSignalLength(), 2)
RESULT

CHECK(void setSpacing(const double spacing))
  Test cwt;
  cwt.setSpacing(0.2);
  
  TEST_REAL_EQUAL(cwt.getSpacing(), 0.2)
RESULT

CHECK(void setWavelet(const std::vector<double>& wavelet))
  vector<double> w(1);
  w[1] = 0.5;
  
  Test cwt;
  cwt.getWavelet() = w;
  
  TEST_EQUAL(cwt.getWavelet() == w, true)
RESULT

CHECK((void transform(RawDataPointConstIterator begin_input, RawDataPointConstIterator end_input, float resolution)))
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



