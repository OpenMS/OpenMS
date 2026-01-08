// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FEATUREFINDER/MultiplexSatelliteCentroided.h>

using namespace OpenMS;

START_TEST(MultiplexSatelliteCentroided, "$Id$")

MultiplexSatelliteCentroided* nullPointer = nullptr;
MultiplexSatelliteCentroided* ptr;

START_SECTION(MultiplexSatelliteCentroided(size_t rt_idx, size_t mz_idx))
  MultiplexSatelliteCentroided satellite(4, 7);
  TEST_EQUAL(satellite.getMZidx(), 7);
  ptr = new MultiplexSatelliteCentroided(4, 7);
  TEST_NOT_EQUAL(ptr, nullPointer);
  delete ptr;
END_SECTION

START_SECTION(size_t getMZidx())
  MultiplexSatelliteCentroided satellite(4, 7);
  TEST_EQUAL(satellite.getMZidx(), 7);
END_SECTION

START_SECTION(size_t getRTidx())
  MultiplexSatelliteCentroided satellite(4, 7);
  TEST_EQUAL(satellite.getRTidx(), 4);
END_SECTION

END_TEST
