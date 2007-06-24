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


CHECK(static Real getInternalToFullAverageWeight())
	TEST_EQUAL(e_ptr->getInternalToFullAverageWeight(), h2o.getAverageWeight())
RESULT

CHECK(static Real getInternalToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getNTerminalToFull())

RESULT

CHECK(static Real getNTerminalToFullAverageWeight())

RESULT

CHECK(static Real getNTerminalToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getCTerminalToFull())

RESULT

CHECK(static Real getCTerminalToFullAverageWeight())

RESULT

CHECK(static Real getCTerminalToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getBIonToFull())

RESULT

CHECK(static Real getBIonToFullAverageWeight())

RESULT

CHECK(static Real getBIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getAIonToFull())

RESULT

CHECK(static Real getAIonToFullAverageWeight())

RESULT

CHECK(static Real getAIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getYIonToFull())

RESULT

CHECK(static Real getYIonToFullAverageWeight())

RESULT

CHECK(static Real getYIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getCIonToFull())

RESULT

CHECK(static Real getCIonToFullAverageWeight())

RESULT

CHECK(static Real getCIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getXIonToFull())

RESULT

CHECK(static Real getXIonToFullAverageWeight())

RESULT

CHECK(static Real getXIonToFullMonoWeight())

RESULT

CHECK(static const EmpiricalFormula& getZIonToFull())

RESULT

CHECK(static Real getZIonToFullAverageWeight())

RESULT

CHECK(static Real getZIonToFullMonoWeight())

RESULT

END_TEST
