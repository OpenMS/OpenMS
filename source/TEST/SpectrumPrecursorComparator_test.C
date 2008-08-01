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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/COMPARISON/SPECTRA/SpectrumPrecursorComparator.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectrumPrecursorComparator, "$Id$")

/////////////////////////////////////////////////////////////

SpectrumPrecursorComparator* e_ptr = 0;
CHECK(SpectrumPrecursorComparator())
	e_ptr = new SpectrumPrecursorComparator;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~SpectrumPrecursorComparator())
	delete e_ptr;
RESULT

e_ptr = new SpectrumPrecursorComparator();

CHECK(SpectrumPrecursorComparator(const SpectrumPrecursorComparator& source))
	SpectrumPrecursorComparator copy(*e_ptr);
	TEST_EQUAL(copy.getName(), e_ptr->getName())
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
RESULT

CHECK(SpectrumPrecursorComparator& operator = (const SpectrumPrecursorComparator& source))
	SpectrumPrecursorComparator copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getName(), e_ptr->getName())
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
RESULT

CHECK(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load("data/Transformers_tests.dta", spec1);

	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load("data/Transformers_tests_2.dta", spec2);

	double score = (*e_ptr)(spec1, spec2);

	TEST_REAL_EQUAL(score, 1.7685)

	score = (*e_ptr)(spec1, spec1);

	TEST_REAL_EQUAL(score, 2)
RESULT

CHECK(double operator () (const PeakSpectrum& a) const)
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load("data/Transformers_tests.dta", spec1);

	TEST_REAL_EQUAL((*e_ptr)(spec1), 2.0)

RESULT

CHECK(static PeakSpectrumCompareFunctor* create())
	PeakSpectrumCompareFunctor* cf = SpectrumPrecursorComparator::create();
	SpectrumPrecursorComparator pre_comp;
	TEST_EQUAL(cf->getName(), pre_comp.getName())
	TEST_EQUAL(cf->getParameters(), pre_comp.getParameters())
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(e_ptr->getProductName(), "SpectrumPrecursorComparator")
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
