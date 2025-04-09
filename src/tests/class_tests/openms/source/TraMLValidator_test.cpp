// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/VALIDATORS/TraMLValidator.h>
///////////////////////////
#include <OpenMS/FORMAT/ControlledVocabulary.h>

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

START_TEST(TraMLValidator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVMappings mapping;
ControlledVocabulary cv;

TraMLValidator* ptr = nullptr;
TraMLValidator* nullPointer = nullptr;
START_SECTION((TraMLValidator(const CVMappings &mapping, const ControlledVocabulary &cv)))
{
	ptr = new TraMLValidator(mapping, cv);
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~TraMLValidator())
{
	delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



