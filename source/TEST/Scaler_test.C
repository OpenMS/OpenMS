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

#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Scaler, "$Id$")

/////////////////////////////////////////////////////////////

PRECISION(0.01)

Scaler* e_ptr = 0;
CHECK((Scaler()))
	e_ptr = new Scaler;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK((~Scaler()))
	delete e_ptr;
RESULT

e_ptr = new Scaler();

CHECK((Scaler(const Scaler& source)))
	Scaler copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((Scaler& operator = (const Scaler& source)))
	Scaler copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/Transformers_tests.dta", spec);

	e_ptr->filterSpectrum(spec);
	
	TEST_EQUAL(spec.size(), 121)

	spec.getContainer().sortByIntensity();
	TEST_REAL_EQUAL(spec.begin()->getIntensity(), 96)
	TEST_REAL_EQUAL((spec.end()-1)->getIntensity(), 121)
	TEST_REAL_EQUAL((spec.end()-1)->getPosition()[0], 136.077)
	
RESULT

CHECK((static PreprocessingFunctor* create()))
	// only tested within the factory
RESULT

CHECK((static const String getName()))
	TEST_EQUAL(e_ptr->getName(), "Scaler")
RESULT

CHECK((void filterPeakMap(PeakMap& exp)))
	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load("data/Transformers_tests.dta", spec);

	PeakMap pm;
	pm.push_back(spec);

  e_ptr->filterPeakMap(pm);

  TEST_EQUAL(pm.begin()->size(), 121)

  pm.begin()->getContainer().sortByIntensity();
  TEST_REAL_EQUAL(pm.begin()->begin()->getIntensity(), 96)
  TEST_REAL_EQUAL((pm.begin()->end()-1)->getIntensity(), 121)
  TEST_REAL_EQUAL((pm.begin()->end()-1)->getPosition()[0], 136.077)

RESULT

CHECK((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load("data/Transformers_tests.dta", spec);

  e_ptr->filterPeakSpectrum(spec);

  TEST_EQUAL(spec.size(), 121)

  spec.getContainer().sortByIntensity();
  TEST_REAL_EQUAL(spec.begin()->getIntensity(), 96)
  TEST_REAL_EQUAL((spec.end()-1)->getIntensity(), 121)
  TEST_REAL_EQUAL((spec.end()-1)->getPosition()[0], 136.077)
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
