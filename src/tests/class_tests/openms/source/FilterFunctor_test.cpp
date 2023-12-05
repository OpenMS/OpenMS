// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

#include <OpenMS/CONCEPT/Factory.h>
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

START_SECTION(static void registerChildren())
	FilterFunctor* ff = Factory<FilterFunctor>::create("ComplementFilter");
	TEST_EQUAL(ff->getName(), "ComplementFilter")
	delete ff;

	ff = Factory<FilterFunctor>::create("IntensityBalanceFilter");
	TEST_EQUAL(ff->getName(), "IntensityBalanceFilter")
	delete ff;

	ff = Factory<FilterFunctor>::create("NeutralLossDiffFilter");
	TEST_EQUAL(ff->getName(), "NeutralLossDiffFilter")
	delete ff;

	ff = Factory<FilterFunctor>::create("IsotopeDiffFilter");
	TEST_EQUAL(ff->getName(), "IsotopeDiffFilter")
	delete ff;

	ff = Factory<FilterFunctor>::create("TICFilter");
	TEST_EQUAL(ff->getName(), "TICFilter")
	delete ff;
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
