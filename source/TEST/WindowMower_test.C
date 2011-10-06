// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(WindowMower, "$Id$")

/////////////////////////////////////////////////////////////

WindowMower* e_ptr = 0;
WindowMower* e_nullPointer = 0;
START_SECTION((WindowMower()))
	e_ptr = new WindowMower;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~WindowMower()))
	delete e_ptr;
END_SECTION

e_ptr = new WindowMower();

START_SECTION((WindowMower(const WindowMower& source)))
	WindowMower copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((WindowMower& operator = (const WindowMower& source)))
	WindowMower copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	TEST_EQUAL(spec.size(), 121)

	Param p(e_ptr->getParameters());
	p.setValue("windowsize", 50.0); // default
	p.setValue("peakcount", 2);
	e_ptr->setParameters(p);
	
	e_ptr->filterSpectrum(spec);
	
	TEST_EQUAL(spec.size(), 56)
	
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.push_back(spec);

  TEST_EQUAL(pm.begin()->size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("windowsize", 50.0); // default
  p.setValue("peakcount", 2);
  e_ptr->setParameters(p);

  e_ptr->filterPeakMap(pm);

  TEST_EQUAL(pm.begin()->size(), 56)
END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
  TEST_EQUAL(spec.size(), 121)

  Param p(e_ptr->getParameters());
  p.setValue("windowsize", 50.0); // default
  p.setValue("peakcount", 2);
  e_ptr->setParameters(p);

  e_ptr->filterPeakSpectrum(spec);

  TEST_EQUAL(spec.size(), 56)
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
