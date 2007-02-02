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

#include <OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ComplementFilter, "$Id$")

/////////////////////////////////////////////////////////////

ComplementFilter* e_ptr = 0;
CHECK((ComplementFilter()))
	e_ptr = new ComplementFilter;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK((~ComplementFilter()))
	delete e_ptr;
RESULT

e_ptr = new ComplementFilter();

CHECK((ComplementFilter(const ComplementFilter& source)))
	ComplementFilter copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((ComplementFilter& operator = (const ComplementFilter& source)))
	ComplementFilter copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((template<typename SpectrumType> double apply(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/Transformers_tests.dta", spec);

	double filter = e_ptr->apply(spec);
	TEST_REAL_EQUAL(filter, 37.0)

	Param p;
	p.setValue("tolerance", 2);
	e_ptr->setParameters(p);
	filter = e_ptr->apply(spec);
	TEST_REAL_EQUAL(filter, 132.5)
	
RESULT

CHECK((static FilterFunctor* create()))
	FilterFunctor* ff = ComplementFilter::create();
	ComplementFilter cf;
	TEST_EQUAL(ff->getParameters(), cf.getParameters())
	TEST_EQUAL(ff->getName(), cf.getName())
RESULT

CHECK((static const String getName()))
	TEST_EQUAL(e_ptr->getName(), "ComplementFilter")
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
