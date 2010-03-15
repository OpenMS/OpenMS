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
// $Maintainer: Chris Bielow, Andreas Bertsch $
// $Authors: Chris Bielow, Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/FILTERING/TRANSFORMERS/SpectraMerger.h>
#include <OpenMS/KERNEL/StandardTypes.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectraMerger, "$Id: SpectraMerger_test.C 6565 2010-01-23 17:21:27Z andreas_bertsch $")

/////////////////////////////////////////////////////////////

SpectraMerger* e_ptr = 0;
START_SECTION((SpectraMerger()))
	e_ptr = new SpectraMerger;
	TEST_NOT_EQUAL(e_ptr, 0)
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

START_SECTION((template <typename ExperimentType> void mergeSpectraBlockWise(ExperimentType& spectrum)))
	// TODO
END_SECTION

START_SECTION((template <typename ExperimentType> void mergeSpectraPrecursors(ExperimentType& spectrum)))
	// TODO
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
