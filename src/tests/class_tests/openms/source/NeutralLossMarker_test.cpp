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

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <map>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(NeutralLossMarker, "$Id$")

/////////////////////////////////////////////////////////////

NeutralLossMarker* e_ptr = nullptr;
NeutralLossMarker* e_nullPointer = nullptr;

START_SECTION((NeutralLossMarker()))
	e_ptr = new NeutralLossMarker;
  TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~NeutralLossMarker()))
	delete e_ptr;
END_SECTION

e_ptr = new NeutralLossMarker();

START_SECTION((NeutralLossMarker(const NeutralLossMarker& source)))
	NeutralLossMarker copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((NeutralLossMarker& operator = (const NeutralLossMarker& source)))
	NeutralLossMarker copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void apply(std::map<double, bool>& marked, SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	map<double, bool> marked;
	e_ptr->apply(marked, spec);

	TEST_EQUAL(marked.size(), 17)
	
	Param p(e_ptr->getParameters());
	p.setValue("tolerance", 10.0);
	e_ptr->setParameters(p);

	marked.clear();
	e_ptr->apply(marked, spec);
	TEST_EQUAL(marked.size(), 49)
END_SECTION

START_SECTION((static PeakMarker* create()))
	PeakMarker* pm = NeutralLossMarker::create();
	NeutralLossMarker marker;
	TEST_EQUAL(pm->getParameters(), marker.getParameters())
	TEST_EQUAL(pm->getName(), marker.getName())
	delete pm;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "NeutralLossMarker")
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
