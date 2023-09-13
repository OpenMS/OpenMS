// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FORMAT/VALIDATORS/MzQuantMLValidator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MzQuantMLValidator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MzQuantMLValidator* ptr = 0;
MzQuantMLValidator* null_ptr = 0;
START_SECTION(MzQuantMLValidator())
{
	ptr = new MzQuantMLValidator();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MzQuantMLValidator())
{
	delete ptr;
}
END_SECTION

START_SECTION((MzQuantMLValidator(const CVMappings &mapping, const ControlledVocabulary &cv)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~MzQuantMLValidator()))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



