// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Volker Mosthaf $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/SqrtMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <cmath>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SqrtMower, "$Id$")

/////////////////////////////////////////////////////////////

SqrtMower* e_ptr = 0;
START_SECTION((SqrtMower()))
	e_ptr = new SqrtMower;
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION((~SqrtMower()))
	delete e_ptr;
END_SECTION

e_ptr = new SqrtMower();

START_SECTION((SqrtMower(const SqrtMower& source)))
	SqrtMower copy(*e_ptr);
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

START_SECTION((SqrtMower& operator=(const SqrtMower& source)))
	SqrtMower copy;
	copy = *e_ptr;
	TEST_EQUAL(*e_ptr == copy, true)
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	
	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), 37.5)

	e_ptr->filterSpectrum(spec);
	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), sqrt(37.5))
END_SECTION

START_SECTION((static PreprocessingFunctor* create()))
	NOT_TESTABLE
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "SqrtMower")
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
	DTAFile dta_file;
  PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.push_back(spec);

	TEST_REAL_SIMILAR((pm.begin()->begin() + 40)->getIntensity(), 37.5)

	e_ptr->filterPeakMap(pm);
	TEST_REAL_SIMILAR((pm.begin()->begin() + 40)->getIntensity(), sqrt(37.5))
END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), 37.5)

	e_ptr->filterPeakSpectrum(spec);
	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), sqrt(37.5))
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
