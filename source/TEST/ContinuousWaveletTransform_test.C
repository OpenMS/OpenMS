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

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransform.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ContinuousWaveletTransform, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ContinuousWaveletTransform* ptr = 0;
ContinuousWaveletTransform* nullPointer = 0;
START_SECTION((ContinuousWaveletTransform()))
  ptr = new ContinuousWaveletTransform();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~ContinuousWaveletTransform()))
  delete ptr;
END_SECTION


START_SECTION((std::vector<Peak1D >& getSignal()))
  ContinuousWaveletTransform cwt;
  
  std::vector<Peak1D > signal(2);
  cwt.getSignal() = signal;
  
  TEST_EQUAL(cwt.getSignal() == signal, true)
END_SECTION

START_SECTION((const std::vector<Peak1D >& getSignal() const))
 ContinuousWaveletTransform cwt;
 
 TEST_EQUAL(cwt.getSignal().empty(), true)
END_SECTION

START_SECTION((DoubleReal getScale() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_SIMILAR(cwt.getScale(), 0)
END_SECTION

START_SECTION((DoubleReal getSpacing() const))
  ContinuousWaveletTransform cwt;
  
  TEST_REAL_SIMILAR(cwt.getSpacing(), 0)
END_SECTION

START_SECTION((double operator[](unsigned int i) const ))
  ContinuousWaveletTransform cwt;
  
  std::vector<Peak1D > signal(1);
  cwt.getSignal() = signal;
  
  ContinuousWaveletTransform const cwt_const(cwt);
  
  TEST_REAL_SIMILAR(cwt_const[0],0)
END_SECTION

START_SECTION((SignedSize getLeftPaddingIndex() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getLeftPaddingIndex(), 0)
END_SECTION

START_SECTION((SignedSize getRightPaddingIndex() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getRightPaddingIndex(), 0)
END_SECTION

START_SECTION((SignedSize getSignalLength() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getSignalLength(), 0)
END_SECTION

START_SECTION((int getSize() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getSize(), 0)
END_SECTION

START_SECTION((const std::vector<double>& getWavelet() const))
  ContinuousWaveletTransform cwt;
  
  TEST_EQUAL(cwt.getWavelet().size(), 0)
END_SECTION

START_SECTION((double& getScale()))
  ContinuousWaveletTransform cwt;
  cwt.getScale() = 0.2;
  
  TEST_REAL_SIMILAR(cwt.getScale(), 0.2)
END_SECTION

START_SECTION((double& getSpacing()))
  ContinuousWaveletTransform cwt;
  cwt.getSpacing() = 0.2;
  
  TEST_REAL_SIMILAR(cwt.getSpacing(), 0.2)
END_SECTION

START_SECTION((double operator[](unsigned int i)))
  std::vector<Peak1D > signal;
  Peak1D rp;
  rp.setIntensity(100.0f);
  signal.push_back(rp);
  
  ContinuousWaveletTransform cwt;
  cwt.getSignal() = signal;
  
  TEST_EQUAL(cwt[0],100)
END_SECTION

START_SECTION((SignedSize& getLeftPaddingIndex()))
  ContinuousWaveletTransform cwt;
  cwt.getLeftPaddingIndex() = 2;
  
  TEST_EQUAL(cwt.getLeftPaddingIndex(), 2)
END_SECTION

START_SECTION((SignedSize& getRightPaddingIndex()))
  ContinuousWaveletTransform cwt;
  cwt.getRightPaddingIndex() = 2;
  
  TEST_EQUAL(cwt.getRightPaddingIndex(), 2)
END_SECTION

START_SECTION((SignedSize& getSignalLength()))
  ContinuousWaveletTransform cwt;
  cwt.getSignalLength() = 2;
  
  TEST_EQUAL(cwt.getSignalLength(), 2)
END_SECTION

START_SECTION((std::vector<double>& getWavelet()))
  vector<double> w(1);
  w[0] = 0.5;
  
  ContinuousWaveletTransform cwt;
  cwt.getWavelet() = w;
  
  TEST_EQUAL(cwt.getWavelet() == w, true)
END_SECTION



START_SECTION((virtual void init(double scale, double spacing)))
  ContinuousWaveletTransform cwt;
  double scale = 0.2;
  double spacing = 2.3;
  cwt.init(scale,spacing);
  
  TEST_REAL_SIMILAR(cwt.getSpacing(),spacing)
  TEST_REAL_SIMILAR(cwt.getScale(),scale)
END_SECTION

START_SECTION((void setLeftPaddingIndex(const SignedSize end_left_padding)))
  ContinuousWaveletTransform cwt;
  cwt.setLeftPaddingIndex(2);
  
  TEST_EQUAL(cwt.getLeftPaddingIndex(), 2)
END_SECTION

START_SECTION((void setRightPaddingIndex(const SignedSize begin_right_padding)))
  ContinuousWaveletTransform cwt;
  cwt.setRightPaddingIndex(2);
  
  TEST_EQUAL(cwt.getRightPaddingIndex(), 2)
END_SECTION

START_SECTION((void setScale(DoubleReal scale)))
  ContinuousWaveletTransform cwt;
  cwt.setScale(0.2);
  
  TEST_REAL_SIMILAR(cwt.getScale(), 0.2)
END_SECTION

START_SECTION((void setSignal(const std::vector<Peak1D >& signal)))
  ContinuousWaveletTransform cwt;
  
  std::vector<Peak1D > signal(2);
  cwt.setSignal(signal);
  
  TEST_EQUAL(cwt.getSignal() == signal, true)
END_SECTION

START_SECTION((void setSignalLength(const SignedSize signal_length)))
  ContinuousWaveletTransform cwt;
  cwt.setSignalLength(2);
  
  TEST_EQUAL(cwt.getSignalLength(), 2)
END_SECTION

START_SECTION((void setSpacing(double spacing)))
  ContinuousWaveletTransform cwt;
  cwt.setSpacing(0.2);
  
  TEST_REAL_SIMILAR(cwt.getSpacing(), 0.2)
END_SECTION

START_SECTION((void setWavelet(const std::vector<double>& wavelet)))
  vector<double> w(1);
  w[0] = 0.5;
  
  ContinuousWaveletTransform cwt;
  cwt.getWavelet() = w;
  
  TEST_EQUAL(cwt.getWavelet() == w, true)
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



