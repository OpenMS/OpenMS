// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow$
// $Authors: Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/TIC.h>

///////////////////////////

START_TEST(TIC, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
using namespace OpenMS;
using namespace std;

TIC* ptr = nullptr;
TIC* nullPointer = nullptr;
START_SECTION(TIC())
ptr = new TIC();
TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~TIC())
delete ptr;
END_SECTION

TIC tic;
START_SECTION(const String& getName() const override)
TEST_EQUAL(tic.getName(), "TIC")
END_SECTION

START_SECTION(Status requirements() const override)
TEST_EQUAL((tic.requirements() == QCBase::Status(QCBase::Requires::RAWMZML)), true);
END_SECTION

START_SECTION(void compute(const MSExperiment& exp, float bin_size))
// very simple test ATM, check if compute returns an empty Result struct
MSExperiment exp;
TEST_EQUAL(tic.compute(exp, 0) == TIC::Result(), true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
