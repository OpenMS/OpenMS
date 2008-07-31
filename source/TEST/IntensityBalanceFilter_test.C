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

#include <OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IntensityBalanceFilter, "$Id$")

/////////////////////////////////////////////////////////////

IntensityBalanceFilter* e_ptr = 0;
CHECK((IntensityBalanceFilter()))
	e_ptr = new IntensityBalanceFilter;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK((~IntensityBalanceFilter()))
	delete e_ptr;
RESULT

e_ptr = new IntensityBalanceFilter();

CHECK((IntensityBalanceFilter(const IntensityBalanceFilter& source)))
	IntensityBalanceFilter copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((IntensityBalanceFilter& operator=(const IntensityBalanceFilter& source)))
	IntensityBalanceFilter copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((template<typename SpectrumType> double apply(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/Transformers_tests.dta", spec);

	double filter = e_ptr->apply(spec);
	TEST_REAL_EQUAL(filter, 0.842697)

RESULT

CHECK((static FilterFunctor* create()))
	FilterFunctor* ff = IntensityBalanceFilter::create();
	IntensityBalanceFilter filter;
	TEST_EQUAL(ff->getParameters(), filter.getParameters())
	TEST_EQUAL(ff->getName(), filter.getName())
RESULT

CHECK((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "IntensityBalanceFilter")
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
