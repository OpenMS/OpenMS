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

#include <OpenMS/CHEMISTRY/Residue.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Residue, "$Id$")

/////////////////////////////////////////////////////////////

Residue* e_ptr = 0;
CHECK(Residue())
	e_ptr = new Residue();
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~Residue())
	delete e_ptr;
RESULT


e_ptr = new Residue();
EmpiricalFormula h2o("H2O");

CHECK(inline static const EmpiricalFormula& getInternalToFull())
	TEST_EQUAL(e_ptr->getInternalToFull(), h2o)
RESULT


CHECK(static DoubleReal getInternalToFullAverageWeight())
	TEST_EQUAL(e_ptr->getInternalToFullAverageWeight(), h2o.getAverageWeight())
RESULT

CHECK(static DoubleReal getInternalToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getNTerminalToFull())

RESULT

CHECK(static DoubleReal getNTerminalToFullAverageWeight())

RESULT

CHECK(static DoubleReal getNTerminalToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getCTerminalToFull())

RESULT

CHECK(static DoubleReal getCTerminalToFullAverageWeight())

RESULT

CHECK(static DoubleReal getCTerminalToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getBIonToFull())

RESULT

CHECK(static DoubleReal getBIonToFullAverageWeight())

RESULT

CHECK(static DoubleReal getBIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getAIonToFull())

RESULT

CHECK(static DoubleReal getAIonToFullAverageWeight())

RESULT

CHECK(static DoubleReal getAIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getYIonToFull())

RESULT

CHECK(static DoubleReal getYIonToFullAverageWeight())

RESULT

CHECK(static DoubleReal getYIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getCIonToFull())

RESULT

CHECK(static DoubleReal getCIonToFullAverageWeight())

RESULT

CHECK(static DoubleReal getCIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getXIonToFull())

RESULT

CHECK(static DoubleReal getXIonToFullAverageWeight())

RESULT

CHECK(static DoubleReal getXIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getZIonToFull())

RESULT

CHECK(static DoubleReal getZIonToFullAverageWeight())

RESULT

CHECK(static DoubleReal getZIonToFullMonoWeight())

RESULT

END_TEST
