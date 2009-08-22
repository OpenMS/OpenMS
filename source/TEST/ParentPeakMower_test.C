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

#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(ParentPeakMower, "$Id$")

/////////////////////////////////////////////////////////////

ParentPeakMower* e_ptr = 0;
START_SECTION((ParentPeakMower()))
	e_ptr = new ParentPeakMower;
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION((~ParentPeakMower()))
	delete e_ptr;
END_SECTION

e_ptr = new ParentPeakMower();

START_SECTION((ParentPeakMower(const ParentPeakMower& source)))
	ParentPeakMower copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((ParentPeakMower& operator = (const ParentPeakMower& source)))
	ParentPeakMower copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION((template<typename SpectrumType> void filterSpectrum(SpectrumType& spectrum)))
	DTAFile dta_file;
	PeakSpectrum spec;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
	spec.setMSLevel(2);
	
	spec.sortByPosition();

	TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), 37.5)

	double window_size(2.0);
	Param p(e_ptr->getParameters());
	p.setValue("window_size", window_size);
	p.setValue("default_charge", 2);
	p.setValue("clean_all_charge_states", (short)1);
	p.setValue("set_to_zero", (short)1);
	e_ptr->setParameters(p);

	e_ptr->filterSpectrum(spec);
	double pre_1_pos(spec.getPrecursors()[0].getMZ() * spec.getPrecursors()[0].getCharge());
	for (Int z = 1; z != spec.getPrecursors()[0].getCharge(); ++z)
	{
		for (PeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
		{	
			if (fabs(it->getPosition()[0] - pre_1_pos / double(z)) <= window_size)
			{
				TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
			}

			// test if NH3 loss is correct removed
			if (fabs(it->getPosition()[0] - (pre_1_pos - 17.0) / double(z)) <= window_size)
			{
				TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
			}

			if (fabs(it->getPosition()[0] - (pre_1_pos - 18.0) / double(z)) <= window_size)
			{
				TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
			}
		}
	}
	
END_SECTION

START_SECTION((static PreprocessingFunctor* create()))
	PreprocessingFunctor* ppf = ParentPeakMower::create();
	ParentPeakMower ppm;
	TEST_EQUAL(ppf->getParameters(), ppm.getParameters())
	TEST_EQUAL(ppf->getName(), ppm.getName())
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(e_ptr->getProductName(), "ParentPeakMower")
END_SECTION

START_SECTION((void filterPeakMap(PeakMap& exp)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);

	PeakMap pm;
	pm.push_back(spec);

  pm.begin()->setMSLevel(2);

  pm.begin()->sortByPosition();

  TEST_REAL_SIMILAR((pm.begin()->begin() + 40)->getIntensity(), 37.5)

  double window_size(2.0);
	Param p(e_ptr->getParameters());
  p.setValue("window_size", window_size);
  p.setValue("default_charge", 2);
  p.setValue("clean_all_charge_states", (short)1);
  p.setValue("set_to_zero", (short)1);
	e_ptr->setParameters(p);

  e_ptr->filterPeakMap(pm);
  double pre_1_pos(pm.begin()->getPrecursors()[0].getMZ() * pm.begin()->getPrecursors()[0].getCharge());
  for (Int z = 1; z != pm.begin()->getPrecursors()[0].getCharge(); ++z)
  {
    for (PeakMap::SpectrumType::ConstIterator it = pm.begin()->begin(); it != pm.begin()->end(); ++it)
    {
      if (fabs(it->getPosition()[0] - pre_1_pos / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }

      // test if NH3 loss is correct removed
      if (fabs(it->getPosition()[0] - (pre_1_pos - 17.0) / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }

      if (fabs(it->getPosition()[0] - (pre_1_pos - 18.0) / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }
    }
  }


END_SECTION

START_SECTION((void filterPeakSpectrum(PeakSpectrum& spectrum)))
  DTAFile dta_file;
  PeakSpectrum spec;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec);
  spec.setMSLevel(2);

  spec.sortByPosition();

  TEST_REAL_SIMILAR((spec.begin() + 40)->getIntensity(), 37.5)

  double window_size(2.0);
	Param p(e_ptr->getParameters());
  p.setValue("window_size", window_size);
  p.setValue("default_charge", 2);
  p.setValue("clean_all_charge_states", (short)1);
  p.setValue("set_to_zero", (short)1);
	e_ptr->setParameters(p);

  e_ptr->filterPeakSpectrum(spec);
  double pre_1_pos(spec.getPrecursors()[0].getMZ() * spec.getPrecursors()[0].getCharge());
  for (Int z = 1; z != spec.getPrecursors()[0].getCharge(); ++z)
  {
    for (PeakSpectrum::ConstIterator it = spec.begin(); it != spec.end(); ++it)
    {
      if (fabs(it->getPosition()[0] - pre_1_pos / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }

      // test if NH3 loss is correct removed
      if (fabs(it->getPosition()[0] - (pre_1_pos - 17.0) / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }

      if (fabs(it->getPosition()[0] - (pre_1_pos - 18.0) / double(z)) <= window_size)
      {
        TEST_REAL_SIMILAR(it->getIntensity(), 0.0);
      }
    }
  }

	
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
