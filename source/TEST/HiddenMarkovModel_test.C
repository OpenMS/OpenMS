// -*- Mode: C++; tab-width: 2; -*-
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
// $Authors: Andreas Bertsch $
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
START_SECTION([EXTRA](HMMState()))
	state_ptr = new HMMState();
	TEST_NOT_EQUAL(state_ptr, 0)
END_SECTION

START_SECTION([EXTRA](virtual ~HMMState()))
	delete state_ptr;
END_SECTION

state_ptr = 0;

START_SECTION([EXTRA](HMMState(const String& name, bool hidden = true)))
	state_ptr = new HMMState("state_name_hidden", true);
	TEST_NOT_EQUAL(state_ptr, 0)
	state_ptr2 = new HMMState("state_name_emitting", false);
	TEST_NOT_EQUAL(state_ptr2, 0)
END_SECTION

START_SECTION([EXTRA](const String& getName() const))
	TEST_EQUAL(state_ptr->getName(), "state_name_hidden");
	TEST_EQUAL(state_ptr2->getName(), "state_name_emitting");
END_SECTION

START_SECTION([EXTRA](bool isHidden() const))
	TEST_EQUAL(state_ptr->isHidden(), true)
	TEST_EQUAL(state_ptr2->isHidden(), false)
END_SECTION

START_SECTION([EXTRA](void setName(const String& name)))
	state_ptr->setName("state_name_hidden2");
	TEST_EQUAL(state_ptr->getName(), "state_name_hidden2")
	state_ptr->setName("state_name_hidden");
END_SECTION

START_SECTION([EXTRA](void setHidden(bool hidden)))
	state_ptr->setHidden(false);
	TEST_EQUAL(state_ptr->isHidden(), false)
	state_ptr->setHidden(true);
	TEST_EQUAL(state_ptr->isHidden(), true)
END_SECTION

START_SECTION([EXTRA](const std::set<HMMState*>& getPredecessorStates() const))
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 0)
END_SECTION

START_SECTION([EXTRA](const std::set<HMMState*>& getSuccessorStates() const))
	TEST_EQUAL(state_ptr->getSuccessorStates().size(), 0)
END_SECTION

START_SECTION([EXTRA](void addPredecessorState(HMMState* state)))
	state_ptr->addPredecessorState(state_ptr2);
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 1);
	TEST_EQUAL(*state_ptr->getPredecessorStates().begin(), state_ptr2);
END_SECTION

START_SECTION([EXTRA](void deletePredecessorState(HMMState* state)))
	state_ptr->deletePredecessorState(state_ptr2);
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 0);
END_SECTION

START_SECTION([EXTRA](void addSuccessorState(HMMState* state)))
  state_ptr->addSuccessorState(state_ptr2);
  TEST_EQUAL(state_ptr->getSuccessorStates().size(), 1);
  TEST_EQUAL(*state_ptr->getSuccessorStates().begin(), state_ptr2);
END_SECTION

START_SECTION([EXTRA](void deleteSuccessorState(HMMState* state)))
  state_ptr->deleteSuccessorState(state_ptr2);
  TEST_EQUAL(state_ptr->getSuccessorStates().size(), 0);
END_SECTION


// Hidden Markov Model Tests
START_SECTION((HiddenMarkovModel()))
	ptr = new HiddenMarkovModel();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~HiddenMarkovModel()))
	delete ptr;
END_SECTION

ptr = new HiddenMarkovModel();

START_SECTION((Size getNumberOfStates() const))
	TEST_EQUAL(ptr->getNumberOfStates(), 0)
END_SECTION

START_SECTION((void addNewState(HMMState* state)))
	ptr->addNewState(state_ptr);
	TEST_EQUAL(ptr->getNumberOfStates(), 1)
	ptr->addNewState(state_ptr2);
	ptr->addNewState(state_ptr3);
END_SECTION

START_SECTION((HMMState* getState(const String& name)))
	TEST_EQUAL(ptr->getState("state_name_hidden"), state_ptr)
END_SECTION

START_SECTION((const HMMState* getState(const String& name) const))
	TEST_EQUAL(ptr->getState("state_name_hidden"), state_ptr)
END_SECTION

START_SECTION((double getTransitionProbability(const String& s1, const String& s2) const))
	TEST_REAL_SIMILAR(ptr->getTransitionProbability("state_name_hidden", "state_name_emitting"), 0.0)
END_SECTION

START_SECTION((double getTransitionProbability(HMMState* s1, HMMState* s2) const))
	TEST_REAL_SIMILAR(ptr->getTransitionProbability(state_ptr, state_ptr2), 0.0)
END_SECTION

START_SECTION((void setTransitionProbability(const String& s1, const String& s2, double prob)))
	ptr->setTransitionProbability("state_name_hidden", "state_name_emitting", 0.3);
	TEST_REAL_SIMILAR(ptr->getTransitionProbability("state_name_hidden", "state_name_emitting"), 0.3)
END_SECTION

START_SECTION((void setTransitionProbability(HMMState* s1, HMMState* s2, double prob)))
	ptr->setTransitionProbability(state_ptr, state_ptr2, 0.4);
	TEST_REAL_SIMILAR(ptr->getTransitionProbability(state_ptr, state_ptr2), 0.4)
END_SECTION

START_SECTION((void addSynonymTransition(const String& name1, const String& name2, const String& synonym1, const String& synonym2)))
	HMMState* s1 = new HMMState("state_name_hidden2");
	HMMState* s2 = new HMMState("state_name_emitting2");
	ptr->addNewState(s1);
	ptr->addNewState(s2);
	ptr->addSynonymTransition("state_name_hidden", "state_name_emitting", "state_name_hidden2", "state_name_emitting2");
END_SECTION

//START_SECTION((void buildSynonyms()))
//	ptr->buildSynonyms();
//	TEST_REAL_SIMILAR(ptr->getTransitionProbability("state_name_hidden2", "state_name_emitting2"), 0.4)
//END_SECTION

START_SECTION((void setInitialTransitionProbability(const String& state, double prob)))
	ptr->setInitialTransitionProbability("state_name_hidden2", 1.0);
END_SECTION

START_SECTION((void setTrainingEmissionProbability(const String& state, double prob)))

END_SECTION

START_SECTION((void setTrainingEmissionProbability(HMMState* state, double prob)))

END_SECTION
			
START_SECTION((void enableTransition(HMMState* s1, HMMState* s2)))

END_SECTION

START_SECTION((void enableTransition(const String& s1, const String& s2)))

END_SECTION

START_SECTION((void disableTransition(HMMState* s1, HMMState* s2)))
				
END_SECTION

START_SECTION((void disableTransition(const String& s1, const String& s2)))
				
END_SECTION

START_SECTION((void disableTransitions()))
	ptr->disableTransitions();
END_SECTION

START_SECTION((void calculateEmissionProbabilities(Map<HMMState*, double>& emission_probs)))

END_SECTION

START_SECTION((void train()))
				
END_SECTION

START_SECTION((void evaluate()))
				
END_SECTION

START_SECTION((void estimateUntrainedTransitions()))

END_SECTION

START_SECTION([EXTRA](HMMState(const HMMState& state)))
	HMMState copy(*state_ptr);
	TEST_EQUAL(copy.getName(), state_ptr->getName())
	TEST_EQUAL(copy.getSuccessorStates().size(), state_ptr->getSuccessorStates().size())
	TEST_EQUAL(copy.getPredecessorStates().size(), state_ptr->getPredecessorStates().size())
	TEST_EQUAL(copy.isHidden(), state_ptr->isHidden())
END_SECTION

START_SECTION((HiddenMarkovModel(const HiddenMarkovModel& hmm_new)))
	HiddenMarkovModel copy(*ptr);
	TEST_EQUAL(copy.getNumberOfStates(), ptr->getNumberOfStates())
END_SECTION

START_SECTION([EXTRA](HMMState& operator = (const HMMState&)))
	HMMState copy;
	copy = *state_ptr;
	TEST_EQUAL(copy.getName(), state_ptr->getName())
	TEST_EQUAL(copy.getSuccessorStates().size(), state_ptr->getSuccessorStates().size())
	TEST_EQUAL(copy.getPredecessorStates().size(), state_ptr->getPredecessorStates().size())
	TEST_EQUAL(copy.isHidden(), state_ptr->isHidden())
END_SECTION

START_SECTION((HiddenMarkovModel& operator = (const HiddenMarkovModel&)))
	HiddenMarkovModel copy;
	copy = *ptr;
	TEST_EQUAL(copy.getNumberOfStates(), ptr->getNumberOfStates())
END_SECTION

START_SECTION((void clearInitialTransitionProbabilities()))
	ptr->clearInitialTransitionProbabilities();
END_SECTION

START_SECTION((void clearTrainingEmissionProbabilities()))
	ptr->clearTrainingEmissionProbabilities();
END_SECTION

START_SECTION((void dump()))
	NOT_TESTABLE
END_SECTION

START_SECTION((void forwardDump()))
	NOT_TESTABLE
END_SECTION

START_SECTION((void write(std::ostream& out) const))
	stringstream ss;
	ptr->write(ss);
END_SECTION

START_SECTION((void writeGraphMLFile(const String& filename)))
	NOT_TESTABLE // just for convenience provided
END_SECTION

START_SECTION((void setVariableModifications(const StringList &modifications)))

END_SECTION

START_SECTION((void clear()))
	ptr->clear();
	TEST_EQUAL(ptr->getNumberOfStates(), 0)
END_SECTION

START_SECTION(void addNewState(const String &name))
	ptr->addNewState("new_fancy_state");
	TEST_EQUAL(ptr->getNumberOfStates(), 1)
END_SECTION
  
START_SECTION(void setPseudoCounts(double pseudo_counts))
	ptr->setPseudoCounts(10e-3);
END_SECTION

START_SECTION(double getPseudoCounts() const)
	TEST_EQUAL(ptr->getPseudoCounts(), 10e-3)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

