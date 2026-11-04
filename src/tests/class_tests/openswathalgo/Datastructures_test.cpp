// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include "OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/DataStructures.h"
#include "OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h"

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
  s.getDataArrays().back()->description = "Ion mOBILITY (blah)"; // wrong case
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, true)

  // BrukerTimsFile via IMDataConverter (CV MS:1003008)
  s.getDataArrays().back()->description = "raw inverse reduced ion mobility array";
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)
  // IMDataConverter for millisecond IM (CV MS:1002816)
  s.getDataArrays().back()->description = "mean ion mobility array";
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)
  // MSConvert legacy
  s.getDataArrays().back()->description = "inverse reduced ion mobility";
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)
  // ProteoWizard diaPASEF
  s.getDataArrays().back()->description = "mean inverse reduced ion mobility array";
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)
  // PeakPickerIM output
  s.getDataArrays().back()->description = "Ion Mobility Centroid";
  TEST_EQUAL (s.getDriftTimeArray() == nullptr, false)
}
END_SECTION

START_SECTION(ISpectrumAccess_filterByDrift)
{
  auto input = std::make_shared<Spectrum>();
  TEST_EXCEPTION(std::invalid_argument, ISpectrumAccess::filterByDrift({}, 1, 2))
  TEST_EXCEPTION(std::invalid_argument, ISpectrumAccess::filterByDrift(input, 1, 2))
  auto mobility = std::make_shared<BinaryDataArray>();
  mobility->description = "raw inverse reduced ion mobility array";
  input->getDataArrays().push_back(mobility);
  TEST_TRUE(ISpectrumAccess::filterByDrift(input, 1, 2)->getMZArray()->data.empty())

  // Keep both interval edges, and do not change the input.
  input->getMZArray()->data = {400, 500, 600};
  input->getIntensityArray()->data = {10, 20, 30};
  mobility->data = {1, 2, 3};
  const auto output = ISpectrumAccess::filterByDrift(input, 1, 2);
  TEST_EQUAL(output->getMZArray()->data.size(), 2)
  TEST_EQUAL(output->getMZArray()->data[0], 400)
  TEST_EQUAL(output->getMZArray()->data[1], 500)
  TEST_EQUAL(output->getIntensityArray()->data[1], 20)
  TEST_EQUAL(output->getDriftTimeArray()->data[1], 2)
  TEST_EQUAL(input->getMZArray()->data.size(), 3)
  output->getMZArray()->data[0] = -1;
  TEST_EQUAL(input->getMZArray()->data[0], 400)

  mobility->data.pop_back();
  TEST_EXCEPTION(std::invalid_argument, ISpectrumAccess::filterByDrift(input, 1, 2))
  mobility->data.push_back(3);
  input->getIntensityArray()->data.clear();
  TEST_EXCEPTION(std::invalid_argument, ISpectrumAccess::filterByDrift(input, 1, 2))
  input->getDataArrays()[1].reset();
  TEST_EXCEPTION(std::invalid_argument, ISpectrumAccess::filterByDrift(input, 1, 2))
  input->getDataArrays().clear();
  TEST_EXCEPTION(std::invalid_argument, ISpectrumAccess::filterByDrift(input, 1, 2))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
