// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/NeutralLossMarker.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

#include <map>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(NeutralLossMarker, "$Id: NeutralLossMarker_test.C 5908 2009-08-26 13:44:26Z marc_sturm $")

/////////////////////////////////////////////////////////////

NeutralLossMarker* e_ptr = 0;
START_SECTION((NeutralLossMarker()))
	e_ptr = new NeutralLossMarker;
	TEST_NOT_EQUAL(e_ptr, 0)
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
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "NeutralLossMarker")
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
