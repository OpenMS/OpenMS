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
#include <OpenMS/FORMAT/VALIDATORS/MzDataValidator.h>
///////////////////////////

#include <OpenMS/FORMAT/ControlledVocabulary.h>

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

START_TEST(MzDataValidator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CVMappings mapping;
ControlledVocabulary cv;

MzDataValidator* ptr = nullptr;
MzDataValidator* nullPointer = nullptr;
START_SECTION((MzDataValidator(const CVMappings &mapping, const ControlledVocabulary &cv)))
{
	ptr = new MzDataValidator(mapping, cv);
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~MzDataValidator())
{
	delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



