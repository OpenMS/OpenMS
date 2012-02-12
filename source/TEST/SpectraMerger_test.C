// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectraMerger, "$Id$")

/////////////////////////////////////////////////////////////

SpectraMerger* e_ptr = 0;
SpectraMerger* e_nullPointer = 0;
START_SECTION((SpectraMerger()))
	e_ptr = new SpectraMerger;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION((~SpectraMerger()))
	delete e_ptr;
END_SECTION

e_ptr = new SpectraMerger();

START_SECTION((SpectraMerger(const SpectraMerger& source)))
	SpectraMerger copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
END_SECTION

START_SECTION((SpectraMerger& operator=(const SpectraMerger& source)))
	SpectraMerger copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
END_SECTION

START_SECTION((template < typename MapType > void mergeSpectraBlockWise(MapType &exp)))
	PeakMap exp, exp2;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_input_2.mzML"), exp);
	TEST_EQUAL(exp.size(), 144)

  exp2 = exp;

	SpectraMerger merger;
  Param p;
  p.setValue("block_method:rt_block_size", 5);
  p.setValue("block_method:ms_levels", IntList::create(StringList::create("1")));
  merger.setParameters(p);
  merger.mergeSpectraBlockWise(exp);
  TEST_EQUAL(exp.size(), 130);
  exp=exp2;

  p.setValue("block_method:rt_block_size", 4);
  p.setValue("block_method:ms_levels", IntList::create(StringList::create("2")));
  merger.setParameters(p);
  merger.mergeSpectraBlockWise(exp);
  TEST_EQUAL(exp.size(), 50);
  TEST_REAL_SIMILAR(exp[0].getRT(),201.0275)
  TEST_REAL_SIMILAR(exp[1].getRT(),204.34075)
  TEST_EQUAL(exp[1].getMSLevel(), 2);
  TEST_EQUAL(exp[2].getMSLevel(), 1);
  exp=exp2;

  p.setValue("block_method:rt_block_size", 4);
  p.setValue("block_method:ms_levels", IntList::create(StringList::create("1,2")));
  merger.setParameters(p);
  merger.mergeSpectraBlockWise(exp);
  TEST_EQUAL(exp.size(), 37);

END_SECTION

START_SECTION((template < typename MapType > void mergeSpectraPrecursors(MapType &exp)))
	PeakMap exp;
	MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_input_precursor.mzML"), exp);

	SpectraMerger merger;
	TEST_EQUAL(exp.size(), 17)

  Param p;
  p.setValue("mz_binning_width", 0.3, "Max m/z distance of two peaks to be merged.", StringList::create("advanced"));

  p.setValue("mz_binning_width_unit", "Da", "Unit in which the distance between two peaks is given.", StringList::create("advanced"));

  // same precursor MS/MS merging
 	p.setValue("precursor_method:mz_tolerance", 10e-5, "Max m/z distance of the precursor entries of two spectra to be merged in [Da].");
  p.setValue("precursor_method:rt_tolerance", 5.0, "Max RT distance of the precursor entries of two spectra to be merged in [s].");
  merger.setParameters(p);
  merger.mergeSpectraPrecursors(exp);

  PeakMap exp2;
  MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("SpectraMerger_output_precursor.mzML"), exp2);

	TEST_EQUAL(exp.size(), exp2.size());
  ABORT_IF(exp.size() != exp2.size());

  for (Size i=0;i<exp.size();++i)
  {
    TEST_EQUAL(exp[i].size(), exp2[i].size())
    TEST_EQUAL(exp[i].getMSLevel (), exp2[i].getMSLevel ())
  }

END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
