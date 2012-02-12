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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ContinuousWaveletTransformNumIntegration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ContinuousWaveletTransformNumIntegration* ptr = 0;
ContinuousWaveletTransformNumIntegration* nullPointer = 0;
START_SECTION((ContinuousWaveletTransformNumIntegration()))
  ptr = new ContinuousWaveletTransformNumIntegration();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~ContinuousWaveletTransformNumIntegration()))
  delete ptr;
END_SECTION

START_SECTION((virtual void init(double scale, double spacing)))
  ContinuousWaveletTransformNumIntegration transformer;
  float scale = 0.5f;
  float spacing = 0.1f;
  
  transformer.init(scale,spacing);
  TEST_REAL_SIMILAR(transformer.getWavelet()[0],1.)
  TEST_REAL_SIMILAR(transformer.getScale(),scale)
  TEST_REAL_SIMILAR(transformer.getSpacing(),spacing)
END_SECTION

START_SECTION((template <typename InputPeakIterator> void transform(InputPeakIterator begin_input, InputPeakIterator end_input, float resolution, unsigned int zeros=0)))
  ContinuousWaveletTransformNumIntegration transformer;
  float scale = 0.5f;
  float spacing = 0.1f;
  
  transformer.init(scale,spacing);
  std::vector<Peak1D > raw_data(9);
  raw_data[4].setIntensity(1.0f);
  transformer.transform(raw_data.begin(),raw_data.end(),1.);
  TEST_REAL_SIMILAR(transformer[4],0)
  TEST_REAL_SIMILAR(transformer.getWavelet()[0],1.)
  TEST_REAL_SIMILAR(transformer.getScale(),scale)
  TEST_REAL_SIMILAR(transformer.getSpacing(),spacing)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



