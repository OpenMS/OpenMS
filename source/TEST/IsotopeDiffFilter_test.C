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

#include <OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(IsotopeDiffFilter, "$Id$")

/////////////////////////////////////////////////////////////

IsotopeDiffFilter* e_ptr = 0;
START_SECTION((IsotopeDiffFilter()))
	e_ptr = new IsotopeDiffFilter;
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION((~IsotopeDiffFilter()))
	delete e_ptr;
END_SECTION

e_ptr = new IsotopeDiffFilter();

START_SECTION((IsotopeDiffFilter(const IsotopeDiffFilter& source)))
	IsotopeDiffFilter copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((IsotopeDiffFilter& operator = (const IsotopeDiffFilter& source)))
	IsotopeDiffFilter copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> double apply(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	double filter = e_ptr->apply(spec);
	TEST_REAL_SIMILAR(filter, 0.0)

	Param p(e_ptr->getParameters());
	p.setValue("tolerance", 10.0);
	e_ptr->setParameters(p);
	filter = e_ptr->apply(spec);
	TEST_REAL_SIMILAR(filter, 2162)
END_SECTION

START_SECTION((static FilterFunctor* create()))
	FilterFunctor* ff = IsotopeDiffFilter::create();
	IsotopeDiffFilter filter;
	TEST_EQUAL(ff->getParameters(), filter.getParameters())
	TEST_EQUAL(ff->getName(), filter.getName())
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "IsotopeDiffFilter")
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
