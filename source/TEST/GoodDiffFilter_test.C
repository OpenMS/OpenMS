// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(GoodDiffFilter, "$Id$")

/////////////////////////////////////////////////////////////

GoodDiffFilter* e_ptr = 0;
CHECK((GoodDiffFilter()))
	e_ptr = new GoodDiffFilter;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK((~GoodDiffFilter()))
	delete e_ptr;
RESULT

e_ptr = new GoodDiffFilter();

CHECK((GoodDiffFilter(const GoodDiffFilter& source)))
	GoodDiffFilter copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((GoodDiffFilter& operator=(const GoodDiffFilter& source)))
	GoodDiffFilter copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((template<typename SpectrumType> double apply(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/Transformers_tests.dta", spec);

	double filter = e_ptr->apply(spec);

	TEST_REAL_EQUAL(filter, 0.104879)

	Param p(e_ptr->getParameters());
	p.setValue("tolerance", 10);
	e_ptr->setParameters(p);
	filter = e_ptr->apply(spec);
	
	TEST_REAL_EQUAL(filter, 0.811684)
RESULT

CHECK((static FilterFunctor* create()))
	FilterFunctor* ff = GoodDiffFilter::create();
	GoodDiffFilter good;
	TEST_EQUAL(ff->getParameters(), good.getParameters())
	TEST_EQUAL(ff->getName(), good.getName())
RESULT

CHECK((static const String getName()))
	TEST_EQUAL(e_ptr->getName(), "GoodDiffFilter")
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
