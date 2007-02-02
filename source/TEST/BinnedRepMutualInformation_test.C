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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/COMPARISON/SPECTRA/BinnedRepMutualInformation.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(BinnedRepMutualInformation, "$Id$")

/////////////////////////////////////////////////////////////

BinnedRepMutualInformation* e_ptr = 0;
CHECK(BinnedRepMutualInformation())
	e_ptr = new BinnedRepMutualInformation;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~BinnedRepMutualInformation())
	delete e_ptr;
RESULT

e_ptr = new BinnedRepMutualInformation();

CHECK(BinnedRepMutualInformation(const BinnedRepMutualInformation& source))
	BinnedRepMutualInformation copy(*e_ptr);
	TEST_EQUAL(copy.getName(), e_ptr->getName())
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
RESULT

CHECK(BinnedRepMutualInformation& operator = (const BinnedRepMutualInformation& source))
	BinnedRepMutualInformation copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getName(), e_ptr->getName())
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
RESULT
			  

CHECK(double operator () (const BinnedRep& a, const BinnedRep& b) const)
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load("data/Transformers_tests.dta", spec1);
	BinnedRep br1(spec1, 1.0, 1);
	
	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load("data/Transformers_tests_2.dta", spec2);
	BinnedRep br2(spec2, 1.0, 1);
	
	double score = (*e_ptr)(br1, br2);
	PRECISION(0.01)
	TEST_REAL_EQUAL(score, 0.0322523)

RESULT

CHECK(double operator () (const BinnedRep& a) const)
	DTAFile dta_file;
  PeakSpectrum spec1;
  dta_file.load("data/Transformers_tests.dta", spec1);
  BinnedRep br1(spec1, 1.0, 1);

	PRECISION(0.1) /// @todo fix the 32/64 bit precision here (Andreas)
	TEST_REAL_EQUAL((*e_ptr)(br1), 0.986786)
RESULT

CHECK(static BinnedRepCompareFunctor* create())
	BinnedRepCompareFunctor* cf = BinnedRepMutualInformation::create();
	BinnedRepMutualInformation mi;
	TEST_EQUAL(cf->getName(), mi.getName())
	TEST_EQUAL(cf->getParameters(), mi.getParameters())
RESULT

CHECK(static const String getProductName())
	TEST_EQUAL(e_ptr->getProductName(), "BinnedRepMutualInformation")
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
