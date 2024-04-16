// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

///////////////////////////

#include <vector>
#include <iostream>

///////////////////////////
START_TEST(FilterFunctor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

// pure interface, hardly testable

START_SECTION(FilterFunctor())
	NOT_TESTABLE
END_SECTION

START_SECTION(FilterFunctor(const FilterFunctor& source))
	NOT_TESTABLE
END_SECTION

START_SECTION(FilterFunctor& operator = (const FilterFunctor& source))
	NOT_TESTABLE
END_SECTION

START_SECTION(template<typename SpectrumType> double apply(SpectrumType&))
	NOT_TESTABLE
END_SECTION

START_SECTION(~FilterFunctor())
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
