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
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Scaler, "$Id$")

/////////////////////////////////////////////////////////////

TOLERANCE_ABSOLUTE(0.01)

Scaler* e_ptr = 0;
START_SECTION((Scaler()))
	e_ptr = new Scaler;
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION((~Scaler()))
	delete e_ptr;
END_SECTION

e_ptr = new Scaler();

START_SECTION((Scaler(const Scaler& source)))
	Scaler copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((Scaler& operator = (const Scaler& source)))
	Scaler copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	e_ptr->filterSpectrum(spec);
	
	TEST_EQUAL(spec.size(), 121)

	spec.sortByIntensity();
	TEST_REAL_SIMILAR(spec.begin()->getIntensity(), 96)
	TEST_REAL_SIMILAR((spec.end()-1)->getIntensity(), 121)
	TEST_REAL_SIMILAR((spec.end()-1)->getPosition()[0], 136.077)
	
END_SECTION

START_SECTION((static PreprocessingFunctor* create()))	
	NOT_TESTABLE
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "Scaler")
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.push_back(spec);

  e_ptr->filterPeakMap(pm);

  TEST_EQUAL(pm.begin()->size(), 121)

  pm.begin()->sortByIntensity();
  TEST_REAL_SIMILAR(pm.begin()->begin()->getIntensity(), 96)
  TEST_REAL_SIMILAR((pm.begin()->end()-1)->getIntensity(), 121)
  TEST_REAL_SIMILAR((pm.begin()->end()-1)->getPosition()[0], 136.077)

END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  e_ptr->filterPeakSpectrum(spec);

  TEST_EQUAL(spec.size(), 121)

  spec.sortByIntensity();
  TEST_REAL_SIMILAR(spec.begin()->getIntensity(), 96)
  TEST_REAL_SIMILAR((spec.end()-1)->getIntensity(), 121)
  TEST_REAL_SIMILAR((spec.end()-1)->getPosition()[0], 136.077)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
