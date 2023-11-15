// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(PeakMarker, "$Id$")

/////////////////////////////////////////////////////////////

PeakMarker* e_ptr = nullptr;
PeakMarker* e_nullPointer = nullptr;

START_SECTION((PeakMarker()))
	e_ptr = new PeakMarker;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~PeakMarker()))
	delete e_ptr;
END_SECTION

e_ptr = new PeakMarker();

START_SECTION((PeakMarker(const PeakMarker& source)))
	PeakMarker copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void apply(std::map<double, bool>&, SpectrumType&)))
	// only the derived classes implement this function properly
	NOT_TESTABLE
END_SECTION

START_SECTION(static const String getProductName())
	TEST_EQUAL(e_ptr->getProductName(), "PeakMarker")
END_SECTION

START_SECTION((PeakMarker& operator = (const PeakMarker& source)))
	PeakMarker copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
  TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
