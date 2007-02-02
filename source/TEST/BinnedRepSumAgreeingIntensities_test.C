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

#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSumAgreeingIntensities.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(BinnedRepSumAgreeingIntensities, "$Id$")

/////////////////////////////////////////////////////////////

BinnedRepSumAgreeingIntensities* e_ptr = 0;
CHECK(BinnedRepSumAgreeingIntensities())
	e_ptr = new BinnedRepSumAgreeingIntensities;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~BinnedRepSumAgreeingIntensities())
	delete e_ptr;
RESULT

e_ptr = new BinnedRepSumAgreeingIntensities();

CHECK(BinnedRepSumAgreeingIntensities(const BinnedRepSumAgreeingIntensities& source))
	BinnedRepSumAgreeingIntensities copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK(BinnedRepSumAgreeingIntensities& operator = (const BinnedRepSumAgreeingIntensities& source))
	BinnedRepSumAgreeingIntensities copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK(double operator () (const BinnedRep& csa, const BinnedRep& csb) const)
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load("data/Transformers_tests.dta", spec1);
	TEST_EQUAL(spec1.size(), 121)


	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load("data/Transformers_tests_2.dta", spec2);

	BinnedRep br1(spec1, 1.0, 1), br2(spec2, 1.0, 1);

	double score = (*e_ptr)(br1, br2);

	TEST_REAL_EQUAL(score, 0.0934876)

	score = (*e_ptr)(br1, br1);

	TEST_REAL_EQUAL(score, 1.0)

RESULT

CHECK(double operator () (const BinnedRep& a) const)
	DTAFile dta_file;
  PeakSpectrum spec1;
  dta_file.load("data/Transformers_tests.dta", spec1);
  TEST_EQUAL(spec1.size(), 121)

	BinnedRep br1(spec1, 1.0, 1);

	double score = (*e_ptr)(br1, br1);

	TEST_REAL_EQUAL(score, 1.0)
RESULT

CHECK(static BinnedRepCompareFunctor* create())
	BinnedRepCompareFunctor* brcf = BinnedRepSumAgreeingIntensities::create();
	BinnedRepSumAgreeingIntensities sai;
	TEST_EQUAL(brcf->getParameters(), sai.getParameters())
	TEST_EQUAL(brcf->getName(), sai.getName())
RESULT

CHECK(static const String getName())
	TEST_EQUAL(e_ptr->getName(), "BinnedRepSumAgreeingIntensities")
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
