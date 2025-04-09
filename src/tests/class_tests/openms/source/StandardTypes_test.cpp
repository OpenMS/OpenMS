// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg  $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/StandardTypes.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

START_TEST(StandardTypes, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION(StandardTypes)

	STATUS("Note: Here we only check whether the typedefs (as such) are fine.")
	NOT_TESTABLE;

END_SECTION

// super duper macro
#define GOOD_TYPEDEF(Type)											\
{																								\
	Type* ptr = 0;																\
  Type* nullPointer = 0;												\
  START_SECTION(Type())																	\
    ptr = new Type;															\
    TEST_NOT_EQUAL(ptr, nullPointer)											\
END_SECTION																					\
  START_SECTION(~Type())																\
    delete ptr;																	\
END_SECTION																					\
}

GOOD_TYPEDEF(PeakSpectrum)
GOOD_TYPEDEF(PeakMap)
GOOD_TYPEDEF(PeakSpectrum)
GOOD_TYPEDEF(PeakMap)

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
