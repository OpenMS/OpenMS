// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include "OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h"

#include <OpenMS/CONCEPT/ClassTest.h>
using namespace OpenMS;
using namespace std;
using namespace OpenSwath;

///////////////////////////

START_TEST(DataStructures, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(OSSpectrum_empty)
{
  OSSpectrum s;

  TEST_EQUAL (s.getMZArray() == nullptr, false)
  TEST_EQUAL (s.getIntensityArray() == nullptr, false)
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, true)

  TEST_EQUAL (s.getMZArray()->data.size(), 0)
  TEST_EQUAL (s.getIntensityArray()->data.size(), 0)
}
END_SECTION

START_SECTION(OSSpectrum_data)
{
  OSSpectrum s;

  BinaryDataArrayPtr mz(new BinaryDataArray);
  mz->data.push_back(1.5);
  BinaryDataArrayPtr inten(new BinaryDataArray);
  inten->data.push_back(100.1);
  BinaryDataArrayPtr im(new BinaryDataArray);
  im->data.push_back(300.1);
  im->description = "Ion Mobility"; // old format

  s.setMZArray(mz);
  s.setIntensityArray(inten);
  s.getDataArrays().push_back(im);

  TEST_EQUAL (s.getMZArray() == nullptr, false)
  TEST_EQUAL (s.getIntensityArray() == nullptr, false)
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)

  TEST_EQUAL (s.getMZArray()->data.size(), 1)
  TEST_EQUAL (s.getIntensityArray()->data.size(), 1)
  TEST_EQUAL (s.getDriftTimeArray()->data.size(), 1)

  TEST_REAL_SIMILAR (s.getMZArray()->data[0], 1.5)
  TEST_REAL_SIMILAR (s.getIntensityArray()->data[0], 100.1)
  TEST_REAL_SIMILAR (s.getDriftTimeArray()->data[0], 300.1)
}
END_SECTION

START_SECTION(OSSpectrum_data_2)
{
  OSSpectrum s;

  BinaryDataArrayPtr mz(new BinaryDataArray);
  mz->data.push_back(1.5);
  BinaryDataArrayPtr inten(new BinaryDataArray);
  inten->data.push_back(100.1);
  BinaryDataArrayPtr im(new BinaryDataArray);
  im->data.push_back(300.1);
  im->description = "Ion Mobility (MS:1002476)"; // new format

  s.setMZArray(mz);
  s.setIntensityArray(inten);
  s.getDataArrays().push_back(im);

  TEST_EQUAL (s.getMZArray() == nullptr, false)
  TEST_EQUAL (s.getIntensityArray() == nullptr, false)
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)

  TEST_EQUAL (s.getMZArray()->data.size(), 1)
  TEST_EQUAL (s.getIntensityArray()->data.size(), 1)
  TEST_EQUAL (s.getDriftTimeArray()->data.size(), 1)

  TEST_REAL_SIMILAR (s.getMZArray()->data[0], 1.5)
  TEST_REAL_SIMILAR (s.getIntensityArray()->data[0], 100.1)
  TEST_REAL_SIMILAR (s.getDriftTimeArray()->data[0], 300.1)

  s.getDataArrays().back()->description = "";
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, true)
  s.getDataArrays().back()->description = "Ion Mobility (blah)";
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)
  s.getDataArrays().back()->description = "Ion mOBILITY (blah)"; // wrong
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
