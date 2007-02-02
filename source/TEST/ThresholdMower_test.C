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

#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ThresholdMower, "$Id$")

/////////////////////////////////////////////////////////////

ThresholdMower* e_ptr = 0;
CHECK((ThresholdMower()))
	e_ptr = new ThresholdMower;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK((~ThresholdMower()))
	delete e_ptr;
RESULT

e_ptr = new ThresholdMower();

CHECK((ThresholdMower(const ThresholdMower& source)))
	ThresholdMower copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
RESULT

CHECK((ThresholdMower& operator=(const ThresholdMower& source)))
	ThresholdMower copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName());
RESULT

CHECK((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load("data/Transformers_tests.dta", spec);
	
	TEST_EQUAL(spec.size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("threshold", 1);
	e_ptr->setParameters(p);

	e_ptr->filterSpectrum(spec);
	TEST_EQUAL(spec.size(), 121)

	p.setValue("threshold", 10);
	e_ptr->setParameters(p);

	e_ptr->filterSpectrum(spec);
	TEST_EQUAL(spec.size(), 14)
RESULT

CHECK((static PreprocessingFunctor* create()))
	PreprocessingFunctor* ppf = ThresholdMower::create();
	ThresholdMower mower;
	TEST_EQUAL(ppf->getParameters(), mower.getParameters())
	TEST_EQUAL(ppf->getName(), mower.getName())
RESULT

CHECK((static const String getName()))
	TEST_EQUAL(e_ptr->getName(), "ThresholdMower")
RESULT

CHECK((void filterPeakMap(PeakMap& exp)))
	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load("data/Transformers_tests.dta", spec);

	PeakMap pm;
	pm.push_back(spec);

  TEST_EQUAL(pm.begin()->size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("threshold", 1);
	e_ptr->setParameters(p);

  e_ptr->filterPeakMap(pm);
  TEST_EQUAL(pm.begin()->size(), 121)

  p.setValue("threshold", 10);
	e_ptr->setParameters(p);
  e_ptr->filterPeakMap(pm);
  TEST_EQUAL(pm.begin()->size(), 14)

RESULT

CHECK((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load("data/Transformers_tests.dta", spec);

  TEST_EQUAL(spec.size(), 121)

	Param p(e_ptr->getParameters());
 	p.setValue("threshold", 1);
	e_ptr->setParameters(p);

  e_ptr->filterPeakSpectrum(spec);
  TEST_EQUAL(spec.size(), 121)

  p.setValue("threshold", 10);
	e_ptr->setParameters(p);
  e_ptr->filterPeakSpectrum(spec);
  TEST_EQUAL(spec.size(), 14)
RESULT

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
