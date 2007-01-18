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

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ContinuousWaveletTransform, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ContinuousWaveletTransform* ptr = 0;
CHECK((double& getScale()))
  ptr = new ContinuousWaveletTransform();
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((double& getScale()))
  delete ptr;
RESULT

CHECK((ContinuousWaveletTransform& operator=(const ContinuousWaveletTransform& cwt)))
  ContinuousWaveletTransform cwt;
  
  cwt.getScale() = 0.2;
  cwt.getSpacing() = 0.01;
  cwt.getLeftPaddingIndex() = 3;
  cwt.getRightPaddingIndex() = 1;
  cwt.getSignalLength() = 2;
  
  ContinuousWaveletTransform cwt_copy;
  cwt_copy = cwt;
  
  TEST_REAL_EQUAL(cwt_copy.getScale(), 0.2)
  TEST_REAL_EQUAL(cwt_copy.getSpacing(), 0.01)
  TEST_REAL_EQUAL(cwt_copy.getLeftPaddingIndex(),3)
  TEST_REAL_EQUAL(cwt_copy.getRightPaddingIndex(),1)
  TEST_REAL_EQUAL(cwt_copy.getSignalLength(), 2)
RESULT

CHECK((void setWavelet(const std::vector<double>& wavelet)))
  ContinuousWaveletTransform cwt;
  
  cwt.getScale() = 0.2;
  cwt.getSpacing() = 0.01;
  cwt.getLeftPaddingIndex() = 3;
  cwt.getRightPaddingIndex() = 1;
  cwt.getSignalLength() = 2;
  
  ContinuousWaveletTransform cwt_copy(cwt);
  
  TEST_REAL_EQUAL(cwt_copy.getScale(), 0.2)
  TEST_REAL_EQUAL(cwt_copy.getSpacing(), 0.01)
  TEST_REAL_EQUAL(cwt_copy.getLeftPaddingIndex(),3)
  TEST_REAL_EQUAL(cwt_copy.getRightPaddingIndex(),1)
  TEST_REAL_EQUAL(cwt_copy.getSignalLength(), 2)
RESULT

CHECK((DPeakArray<1, DRawDataPoint<1>& getSignal()))
  ContinuousWaveletTransform cwt;
  
  DPeakArray<1, DRawDataPoint<1> > signal(2);
  cwt.getSignal() = signal;
  
  TEST_EQUAL(cwt.getSignal() == signal, true)
RESULT

CHECK((const DPeakArray<1, DRawDataPoint<1>& getSignal() const))
 ContinuousWaveletTransform cwt;
 
 TEST_EQUAL(cwt.getSignal().size() == 0, true)
RESULT

CHECK((const double& getScale() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_EQUAL(cwt.getScale(), 0)
RESULT

CHECK((const double& getSpacing() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_EQUAL(cwt.getSpacing(), 0)
RESULT

CHECK((const double& operator [] (const unsigned int i) const))
  ContinuousWaveletTransform cwt;
  
  DPeakArray<1, DRawDataPoint<1> > signal(1);
  cwt.getSignal() = signal;
  
  ContinuousWaveletTransform const cwt_const(cwt);
  
  TEST_REAL_EQUAL(cwt_const[0],0)
RESULT

CHECK((const int& getLeftPaddingIndex() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_EQUAL(cwt.getLeftPaddingIndex(), 0)
RESULT

CHECK((const int& getRightPaddingIndex() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_EQUAL(cwt.getRightPaddingIndex(), 0)
RESULT

CHECK((const int& getSignalLength() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_EQUAL(cwt.getSignalLength(), 0)
RESULT

CHECK((int getSize() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_EQUAL(cwt.getSize(), 0)
RESULT

CHECK((const std::vector<double>& getWavelet() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_EQUAL(cwt.getWavelet().size(), 0)
RESULT

CHECK((const double& operator [] (const unsigned int i) const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_EQUAL(cwt.getScale(), 0)
RESULT

CHECK((double& getScale()))
  ContinuousWaveletTransform cwt;
  cwt.getScale() = 0.2;
  
  TEST_REAL_EQUAL(cwt.getScale(), 0.2)
RESULT

CHECK((double& getSpacing()))
  ContinuousWaveletTransform cwt;
  cwt.getSpacing() = 0.2;
  
  TEST_REAL_EQUAL(cwt.getSpacing(), 0.2)
RESULT

CHECK((double& operator [] (const unsigned int i)))
  DPeakArray<1, DRawDataPoint<1> > signal;
  DRawDataPoint<1> rp;
  rp.setIntensity(100);
  signal.push_back(rp);
  
  ContinuousWaveletTransform cwt;
  cwt.getSignal() = signal;
  
  TEST_EQUAL(cwt[0],100)
RESULT

CHECK((int& getLeftPaddingIndex()))
  ContinuousWaveletTransform cwt;
  cwt.getLeftPaddingIndex() = 2;
  
  TEST_REAL_EQUAL(cwt.getLeftPaddingIndex(), 2)
RESULT

CHECK((int& getRightPaddingIndex()))
  ContinuousWaveletTransform cwt;
  cwt.getRightPaddingIndex() = 2;
  
  TEST_REAL_EQUAL(cwt.getRightPaddingIndex(), 2)
RESULT

CHECK((int& getSignalLength()))
  ContinuousWaveletTransform cwt;
  cwt.getSignalLength() = 2;
  
  TEST_REAL_EQUAL(cwt.getSignalLength(), 2)
RESULT

CHECK((std::vector<double>& getWavelet()))
  vector<double> w(1);
  w[0] = 0.5;
  
  ContinuousWaveletTransform cwt;
  cwt.getWavelet() = w;
  
  TEST_EQUAL(cwt.getWavelet() == w, true)
RESULT

CHECK((double& operator [] (const unsigned int i)))
  DPeakArray<1, DRawDataPoint<1> > signal;
  signal.resize(1);
  signal[0].getIntensity() = 1;
  
  ContinuousWaveletTransform cwt;
  cwt.setSignal(signal);
  
  TEST_REAL_EQUAL(cwt[0], 1)
RESULT

CHECK((void init(double scale, double spacing)))
  ContinuousWaveletTransform cwt;
  double scale = 0.2;
  double spacing = 2.3;
  cwt.init(scale,spacing);
  
  TEST_REAL_EQUAL(cwt.getSpacing(),spacing)
  TEST_REAL_EQUAL(cwt.getScale(),scale)
RESULT

CHECK((void setLeftPaddingIndex(const int end_left_padding)))
  ContinuousWaveletTransform cwt;
  cwt.setLeftPaddingIndex(2);
  
  TEST_REAL_EQUAL(cwt.getLeftPaddingIndex(), 2)
RESULT

CHECK((void setRightPaddingIndex(const int begin_right_padding)))
  ContinuousWaveletTransform cwt;
  cwt.setRightPaddingIndex(2);
  
  TEST_EQUAL(cwt.getRightPaddingIndex(), 2)
RESULT

CHECK((void setScale(const double& scale)))
  ContinuousWaveletTransform cwt;
  cwt.setScale(0.2);
  
  TEST_REAL_EQUAL(cwt.getScale(), 0.2)
RESULT

CHECK((void setSignal(const DPeakArray<1, DRawDataPoint<1> >& signal)))
  ContinuousWaveletTransform cwt;
  
  DPeakArray<1, DRawDataPoint<1> > signal(2);
  cwt.setSignal(signal);
  
  TEST_EQUAL(cwt.getSignal() == signal, true)
RESULT

CHECK((void setSignalLength(const int signal_length)))
  ContinuousWaveletTransform cwt;
  cwt.setSignalLength(2);
  
  TEST_EQUAL(cwt.getSignalLength(), 2)
RESULT

CHECK((void setSpacing(const double spacing)))
  ContinuousWaveletTransform cwt;
  cwt.setSpacing(0.2);
  
  TEST_REAL_EQUAL(cwt.getSpacing(), 0.2)
RESULT

CHECK((void setWavelet(const std::vector<double>& wavelet)))
  vector<double> w(1);
  w[0] = 0.5;
  
  ContinuousWaveletTransform cwt;
  cwt.getWavelet() = w;
  
  TEST_EQUAL(cwt.getWavelet() == w, true)
RESULT

CHECK((void setSignal(const DPeakArray<1, DRawDataPoint<1> >& signal)))
  ContinuousWaveletTransform cwt;
  
  DPeakArray<1, DRawDataPoint<1> > signal(2);
  cwt.setSignal(signal);
  
  TEST_EQUAL(cwt.getSignal() == signal, true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



