// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/VALIDATORS/MzMLValidator.h>
#include <OpenMS/DATASTRUCTURES/CVMappings.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

///////////////////////////

START_TEST(MzMLValidator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Internal;
using namespace std;

CVMappings mapping;
ControlledVocabulary cv;

SemanticValidator* ptr = nullptr;
SemanticValidator* nullPointer = nullptr;
START_SECTION((MzMLValidator(const CVMappings& mapping, const ControlledVocabulary& cv)))
	ptr = new MzMLValidator(mapping,cv);
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MzMLValidator()))
	delete ptr;
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



