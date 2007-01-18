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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <iostream>

#include <OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(HiddenMarkovModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

HiddenMarkovModel* ptr = 0;
HMMState* state_ptr = 0;
HMMState* state_ptr2 = 0;
HMMState* state_ptr3 = new HMMState("base", true);

// Hidden Markov Model State Tests
CHECK(HMMState())
	state_ptr = new HMMState();
	TEST_NOT_EQUAL(state_ptr, 0)
RESULT

CHECK(~HMMState())
	delete state_ptr;
RESULT

state_ptr = 0;

CHECK(HMMState(const String& name, bool hidden = true))
	state_ptr = new HMMState("state_name_hidden", true);
	TEST_NOT_EQUAL(state_ptr, 0)
	state_ptr2 = new HMMState("state_name_emitting", false);
	TEST_NOT_EQUAL(state_ptr2, 0)
RESULT

CHECK(const String& getName() const)
	TEST_EQUAL(state_ptr->getName(), "state_name_hidden");
	TEST_EQUAL(state_ptr2->getName(), "state_name_emitting");
RESULT

CHECK(bool isHidden() const)
	TEST_EQUAL(state_ptr->isHidden(), true)
	TEST_EQUAL(state_ptr2->isHidden(), false)
RESULT

CHECK(void setName(const String& name))
	state_ptr->setName("state_name_hidden2");
	TEST_EQUAL(state_ptr->getName(), "state_name_hidden2")
	state_ptr->setName("state_name_hidden");
RESULT

CHECK(void setHidden(bool hidden))
	state_ptr->setHidden(false);
	TEST_EQUAL(state_ptr->isHidden(), false)
	state_ptr->setHidden(true);
	TEST_EQUAL(state_ptr->isHidden(), true)
RESULT

CHECK(const std::set<HMMState*>& getPredecessorStates() const)
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 0)
RESULT

CHECK(const std::set<HMMState*>& getSuccessorStates() const)
	TEST_EQUAL(state_ptr->getSuccessorStates().size(), 0)
RESULT

CHECK(void addPredecessorState(HMMState* state))
	state_ptr->addPredecessorState(state_ptr2);
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 1);
	TEST_EQUAL(*state_ptr->getPredecessorStates().begin(), state_ptr2);
RESULT

CHECK(void deletePredecessorState(HMMState* state))
	state_ptr->deletePredecessorState(state_ptr2);
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 0);
RESULT

CHECK(void addSuccessorState(HMMState* state))
  state_ptr->addSuccessorState(state_ptr2);
  TEST_EQUAL(state_ptr->getSuccessorStates().size(), 1);
  TEST_EQUAL(*state_ptr->getSuccessorStates().begin(), state_ptr2);
RESULT

CHECK(void deleteSuccessorState(HMMState* state))
  state_ptr->deleteSuccessorState(state_ptr2);
  TEST_EQUAL(state_ptr->getSuccessorStates().size(), 0);
RESULT


// Hidden Markov Model Tests
CHECK(HiddenMarkovModel())
	ptr = new HiddenMarkovModel();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~HiddenMarkovModel())
	delete ptr;
RESULT

ptr = new HiddenMarkovModel();

CHECK(Size getNumberOfStates() const)
	TEST_EQUAL(ptr->getNumberOfStates(), 0)
RESULT

CHECK(void addNewState(HMMState* state))
	ptr->addNewState(state_ptr);
	TEST_EQUAL(ptr->getNumberOfStates(), 1)
	ptr->addNewState(state_ptr2);
	ptr->addNewState(state_ptr3);
RESULT

CHECK(HMMState* getState(const String& name))
	TEST_EQUAL(ptr->getState("state_name_hidden"), state_ptr)
RESULT

CHECK(const HMMState* getState(const String& name) const)
	TEST_EQUAL(ptr->getState("state_name_hidden"), state_ptr)
RESULT

CHECK(double getTransitionProbability(const String& s1, const String& s2) const)
	TEST_REAL_EQUAL(ptr->getTransitionProbability("state_name_hidden", "state_name_emitting"), 0.0)
RESULT

CHECK(double getTransitionProbability(HMMState* s1, HMMState* s2) const)
	TEST_REAL_EQUAL(ptr->getTransitionProbability(state_ptr, state_ptr2), 0.0)
RESULT

CHECK(void setTransitionProbability(const String& s1, const String& s2, double prob))
	ptr->setTransitionProbability("state_name_hidden", "state_name_emitting", 0.3);
	TEST_REAL_EQUAL(ptr->getTransitionProbability("state_name_hidden", "state_name_emitting"), 0.3)
RESULT

CHECK(void setTransitionProbability(HMMState* s1, HMMState* s2, double prob))
	ptr->setTransitionProbability(state_ptr, state_ptr2, 0.4);
	TEST_REAL_EQUAL(ptr->getTransitionProbability(state_ptr, state_ptr2), 0.4)
RESULT

CHECK(void addSynonymTransition(const String& name1, const String& name2, const String& synonym1, const String& synonym2))
	HMMState* s1 = new HMMState("state_name_hidden2");
	HMMState* s2 = new HMMState("state_name_emitting2");
	ptr->addNewState(s1);
	ptr->addNewState(s2);
	ptr->addSynonymTransition("state_name_hidden", "state_name_emitting", "state_name_hidden2", "state_name_emitting2");
RESULT

CHECK(void buildSynonyms())
	ptr->buildSynonyms();
	TEST_REAL_EQUAL(ptr->getTransitionProbability("state_name_hidden2", "state_name_emitting2"), 0.4)
RESULT

CHECK(void setInitialTransitionProbability(const String& state, double prob))

RESULT

CHECK(void setInitialTransitionProbability(HMMState* state, double prob))

RESULT

CHECK(void setTrainingEmissionProbability(const String& state, double prob))

RESULT

CHECK(void setTrainingEmissionProbability(HMMState* state, double prob))

RESULT
			
CHECK(void enableTransition(HMMState* s1, HMMState* s2))

RESULT

CHECK(void enableTransition(const String& s1, const String& s2))

RESULT

CHECK(void disableTransition(HMMState* s1, HMMState* s2))
				
RESULT

CHECK(void disableTransition(const String& s1, const String& s2))
				
RESULT

CHECK(void disableTransitions())
	ptr->disableTransitions();
RESULT

CHECK(void calculateEmissionProbabilities(HashMap<HMMState*, double>& emission_probs))

RESULT

CHECK(void train())
				
RESULT

CHECK(void evaluate())
				
RESULT

CHECK(void estimateUntrainedTransitions())

RESULT

CHECK(HMMState(const HMMState& state))
	HMMState copy(*state_ptr);
	TEST_EQUAL(copy.getName(), state_ptr->getName())
	TEST_EQUAL(copy.getSuccessorStates().size(), state_ptr->getSuccessorStates().size())
	TEST_EQUAL(copy.getPredecessorStates().size(), state_ptr->getPredecessorStates().size())
	TEST_EQUAL(copy.isHidden(), state_ptr->isHidden())
RESULT

CHECK(HiddenMarkovModel(const HiddenMarkovModel& hmm_new))
	HiddenMarkovModel copy(*ptr);
	TEST_EQUAL(copy.getNumberOfStates(), ptr->getNumberOfStates())
RESULT

CHECK(HMMState& operator = (const HMMState&))
	HMMState copy;
	copy = *state_ptr;
	TEST_EQUAL(copy.getName(), state_ptr->getName())
	TEST_EQUAL(copy.getSuccessorStates().size(), state_ptr->getSuccessorStates().size())
	TEST_EQUAL(copy.getPredecessorStates().size(), state_ptr->getPredecessorStates().size())
	TEST_EQUAL(copy.isHidden(), state_ptr->isHidden())
RESULT

CHECK(HiddenMarkovModel& operator = (const HiddenMarkovModel&))
	HiddenMarkovModel copy;
	copy = *ptr;
	TEST_EQUAL(copy.getNumberOfStates(), ptr->getNumberOfStates())
RESULT

CHECK(void clearInitialTransitionProbabilities())
	// TODO
	ptr->clearInitialTransitionProbabilities();
RESULT

CHECK(void clearTrainingEmissionProbabilities())
	// TODO
	ptr->clearTrainingEmissionProbabilities();
RESULT

CHECK(void dump())
	// nothing to test
RESULT

CHECK(void forwardDump())
	// nothing to test
RESULT

CHECK(void write(std::ostream& out))
	// TODO
	stringstream ss;
	ptr->write(ss);

RESULT

CHECK(void writetoYGFFile(const String& filename))
	// TODO
RESULT

CHECK(void clear())
	ptr->clear();
	TEST_EQUAL(ptr->getNumberOfStates(), 0)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
