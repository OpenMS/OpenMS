// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/ANALYSIS/ID/HiddenMarkovModelLight.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(HiddenMarkovModelLight, "$Id:$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

HiddenMarkovModelLight* ptr = 0;
HMMStateLight* state_ptr = 0;
HMMStateLight* state_ptr2 = 0;

// Hidden Markov Model State Tests
CHECK(HMMStateLight())
	state_ptr = new HMMStateLight();
	TEST_NOT_EQUAL(state_ptr, 0)
RESULT

CHECK(~HMMStateLight())
	delete state_ptr;
RESULT

state_ptr = 0;

CHECK(HMMStateLight(const String& name, bool hidden = true))
	state_ptr = new HMMStateLight(27, true);
	TEST_NOT_EQUAL(state_ptr, 0)
	state_ptr2 = new HMMStateLight(123, false);
	TEST_NOT_EQUAL(state_ptr2, 0)
RESULT

CHECK(Size getIdentifier() const)
	TEST_EQUAL(state_ptr->getIdentifier(), 27);
	TEST_EQUAL(state_ptr2->getIdentifier(), 123);
RESULT

CHECK(bool isHidden() const)
	TEST_EQUAL(state_ptr->isHidden(), true)
	TEST_EQUAL(state_ptr2->isHidden(), false)
RESULT

CHECK(void setIdentifier(Size id))
	state_ptr->setIdentifier(1234);
	TEST_EQUAL(state_ptr->getIdentifier(), 1234)
RESULT

CHECK(void setHidden(bool hidden))
	state_ptr->setHidden(false);
	TEST_EQUAL(state_ptr->isHidden(), false)
	state_ptr->setHidden(true);
	TEST_EQUAL(state_ptr->isHidden(), true)
RESULT

CHECK(const std::set<HMMStateLight*>& getPredecessorStates() const)
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 0)
RESULT

CHECK(const std::set<HMMStateLight*>& getSuccessorStates() const)
	TEST_EQUAL(state_ptr->getSuccessorStates().size(), 0)
RESULT

CHECK(const addPredecessorState(HMMStateLight* state))
	state_ptr->addPredecessorState(state_ptr2);
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 1);
	TEST_EQUAL(*state_ptr->getPredecessorStates().begin(), state_ptr2);
RESULT

CHECK(const deletePredecessorState(HMMStateLight* state))
	state_ptr->deletePredecessorState(state_ptr2);
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 0);
RESULT

CHECK(const addSuccessorState(HMMStateLight* state))
  state_ptr->addSuccessorState(state_ptr2);
  TEST_EQUAL(state_ptr->getSuccessorStates().size(), 1);
  TEST_EQUAL(*state_ptr->getSuccessorStates().begin(), state_ptr2);
RESULT

CHECK(const deleteSuccessorState(HMMStateLight* state))
  state_ptr->deleteSuccessorState(state_ptr2);
  TEST_EQUAL(state_ptr->getSuccessorStates().size(), 0);
RESULT


// Hidden Markov Model Tests
CHECK(HiddenMarkovModelLight())
	ptr = new HiddenMarkovModelLight();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(virtual ~HiddenMarkovModelLight())
	delete ptr;
RESULT

ptr = new HiddenMarkovModelLight();



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
