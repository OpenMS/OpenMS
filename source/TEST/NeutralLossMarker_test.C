// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

START_TEST(NeutralLossMarker, "$Id$")

/////////////////////////////////////////////////////////////

NeutralLossMarker* e_ptr = 0;
CHECK((NeutralLossMarker()))
	e_ptr = new NeutralLossMarker;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK((~NeutralLossMarker()))
	delete e_ptr;
RESULT

e_ptr = new NeutralLossMarker();

CHECK((NeutralLossMarker(const NeutralLossMarker& source)))
	NeutralLossMarker copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((NeutralLossMarker& operator = (const NeutralLossMarker& source)))
	NeutralLossMarker copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((template<typename SpectrumType> void apply(std::map<double, bool>& marked, SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/Transformers_tests.dta", spec);

	map<double, bool> marked;
	e_ptr->apply(marked, spec);

	TEST_EQUAL(marked.size(), 17)
	
	Param p(e_ptr->getParameters());
	p.setValue("tolerance", 10.0);
	e_ptr->setParameters(p);

	marked.clear();
	e_ptr->apply(marked, spec);
	TEST_EQUAL(marked.size(), 49)
RESULT

CHECK((static PeakMarker* create()))
	PeakMarker* pm = NeutralLossMarker::create();
	NeutralLossMarker marker;
	TEST_EQUAL(pm->getParameters(), marker.getParameters())
	TEST_EQUAL(pm->getName(), marker.getName())
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "NeutralLossMarker")
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
