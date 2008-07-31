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

#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeMarker.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

#include <map>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IsotopeMarker, "$Id$")

/////////////////////////////////////////////////////////////

IsotopeMarker* e_ptr = 0;
CHECK((IsotopeMarker()))
	e_ptr = new IsotopeMarker;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK((~IsotopeMarker()))
	delete e_ptr;
RESULT

e_ptr = new IsotopeMarker();

CHECK((IsotopeMarker(const IsotopeMarker& source)))
	IsotopeMarker copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((IsotopeMarker& operator=(const IsotopeMarker& source)))
	IsotopeMarker copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((template<typename SpectrumType> void apply(std::map<double, bool> marked, SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/Transformers_tests.dta", spec);

	map<double, bool> marked;
	e_ptr->apply(marked, spec);
	
	TEST_EQUAL(marked.size(), 0)

	/// @todo Fix this (Andreas)
	//e_ptr->getParameters().setValue("n", 10);
	//e_ptr->apply(spec);
	//TEST_EQUAL(spec.size(), 10)
RESULT

CHECK((static PeakMarker* create()))
	PeakMarker* pm = IsotopeMarker::create();
	IsotopeMarker im;
	TEST_EQUAL(pm->getParameters(), im.getParameters())
	TEST_EQUAL(pm->getName(), im.getName())
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "IsotopeMarker")
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
