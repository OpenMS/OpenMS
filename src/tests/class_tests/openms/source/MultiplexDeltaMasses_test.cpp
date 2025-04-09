// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FEATUREFINDER/MultiplexDeltaMasses.h>

using namespace OpenMS;

START_TEST(MultiplexDeltaMasses, "$Id$")

MultiplexDeltaMasses* nullPointer = nullptr;
MultiplexDeltaMasses* ptr;

START_SECTION(MultiplexDeltaMasses())
    MultiplexDeltaMasses pattern;
    TEST_EQUAL(pattern.getDeltaMasses().size(), 0);
    ptr = new MultiplexDeltaMasses();
    TEST_NOT_EQUAL(ptr, nullPointer);
    delete ptr;
END_SECTION

MultiplexDeltaMasses pattern;
pattern.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(0, "no_label"));
pattern.getDeltaMasses().push_back(MultiplexDeltaMasses::DeltaMass(6.031817, "Arg6"));

START_SECTION(std::vector<DeltaMass>& getDeltaMasses())
  TEST_REAL_SIMILAR(pattern.getDeltaMasses()[0].delta_mass, 0);
  TEST_REAL_SIMILAR(pattern.getDeltaMasses()[1].delta_mass, 6.031817);
END_SECTION

END_TEST
