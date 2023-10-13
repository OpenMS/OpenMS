// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/ContinuousWaveletTransformNumIntegration.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ContinuousWaveletTransformNumIntegration, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ContinuousWaveletTransformNumIntegration* ptr = nullptr;
ContinuousWaveletTransformNumIntegration* nullPointer = nullptr;
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



