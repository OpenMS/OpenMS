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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Normalizer, "$Id$")

/////////////////////////////////////////////////////////////

Normalizer* e_ptr = 0;
START_SECTION((Normalizer()))
	e_ptr = new Normalizer;
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION((~Normalizer()))
	delete e_ptr;
END_SECTION

e_ptr = new Normalizer();

START_SECTION((Normalizer(const Normalizer& source)))
	Normalizer copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((Normalizer& operator = (const Normalizer& source)))
	Normalizer copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	spec.sortByIntensity();

	TEST_EQUAL(spec.rbegin()->getIntensity(), 46)

	e_ptr->filterSpectrum(spec);

	spec.sortByIntensity();
	
	TEST_EQUAL(spec.rbegin()->getIntensity(), 1)

	Param p(e_ptr->getParameters());
	p.setValue("method", "to_TIC");
	e_ptr->setParameters(p);
	e_ptr->filterSpectrum(spec);

	double sum(0);
	for (PeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
	{
		sum += it->getIntensity();
	}

	TEST_REAL_SIMILAR(sum, 1.0);	
END_SECTION

START_SECTION((static PreprocessingFunctor* create()))
	PreprocessingFunctor* ppf = Normalizer::create();
	Normalizer norm;
	TEST_EQUAL(ppf->getParameters(), norm.getParameters())
	TEST_EQUAL(ppf->getName(), norm.getName())
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "Normalizer")
END_SECTION
	
START_SECTION((void filterPeakMap(PeakMap& exp)))
	delete e_ptr;
	e_ptr = new Normalizer();

	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.push_back(spec);

  pm.begin()->sortByIntensity();

  TEST_EQUAL(pm.begin()->rbegin()->getIntensity(), 46)

  e_ptr->filterPeakMap(pm);

  pm.begin()->sortByIntensity();

  TEST_EQUAL(pm.begin()->rbegin()->getIntensity(), 1)

	Param p(e_ptr->getParameters());
	p.setValue("method", "to_TIC");
	e_ptr->setParameters(p);
  e_ptr->filterPeakMap(pm);

  double sum(0);
  for (PeakMap::SpectrumType::ConstIterator it = pm.begin()->begin(); it != pm.begin()->end(); ++it)
  {
    sum += it->getIntensity();
  }

  TEST_REAL_SIMILAR(sum, 1.0);	
END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
	delete e_ptr;
	e_ptr = new Normalizer();

	DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

  spec.sortByIntensity();

  TEST_EQUAL(spec.rbegin()->getIntensity(), 46)

  e_ptr->filterPeakSpectrum(spec);

  spec.sortByIntensity();

  TEST_EQUAL(spec.rbegin()->getIntensity(), 1)

	Param p(e_ptr->getParameters());
	p.setValue("method", "to_TIC");
	e_ptr->setParameters(p);
  e_ptr->filterPeakSpectrum(spec);

  double sum(0);
  for (PeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
  {
    sum += it->getIntensity();
  }

  TEST_REAL_SIMILAR(sum, 1.0);	
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
