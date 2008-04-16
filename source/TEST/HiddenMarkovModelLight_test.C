// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

START_TEST(HiddenMarkovModelLight, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

HiddenMarkovModelLight* ptr = 0;
HMMStateLight* state_ptr = 0;
HMMStateLight* state_ptr2 = 0;
HMMStateLight* state_ptr3 = new HMMStateLight(0, true);

// Hidden Markov Model State Tests
CHECK([EXTRA]HMMStateLight())
	state_ptr = new HMMStateLight();
	TEST_NOT_EQUAL(state_ptr, 0)
RESULT

CHECK([EXTRA]HMMStateLight(const HMMStateLight& state))
	HMMStateLight copy(*state_ptr);
	TEST_EQUAL(copy.getIdentifier(), state_ptr->getIdentifier())
	TEST_EQUAL(copy.getSuccessorStates().size(), state_ptr->getSuccessorStates().size())
	TEST_EQUAL(copy.getPredecessorStates().size(), state_ptr->getPredecessorStates().size())
	TEST_EQUAL(copy.isHidden(), state_ptr->isHidden())
RESULT

CHECK([EXTRA]HMMStateLight& operator = (const HMMStateLight&))
	HMMStateLight copy;
	copy = *state_ptr;
	TEST_EQUAL(copy.getIdentifier(), state_ptr->getIdentifier())
	TEST_EQUAL(copy.getSuccessorStates().size(), state_ptr->getSuccessorStates().size())
	TEST_EQUAL(copy.getPredecessorStates().size(), state_ptr->getPredecessorStates().size())
	TEST_EQUAL(copy.isHidden(), state_ptr->isHidden())
RESULT

CHECK([EXTRA] virtual ~HMMStateLight())
  delete state_ptr;
RESULT

state_ptr = 0;


CHECK([EXTRA]HMMStateLight(UInt identifier, bool hidden = true))
	state_ptr = new HMMStateLight(27, true);
	TEST_NOT_EQUAL(state_ptr, 0)
	state_ptr2 = new HMMStateLight(123, false);
	TEST_NOT_EQUAL(state_ptr2, 0)
RESULT

CHECK([EXTRA]UInt getIdentifier() const)
	TEST_EQUAL(state_ptr->getIdentifier(), 27);
	TEST_EQUAL(state_ptr2->getIdentifier(), 123);
RESULT

CHECK([EXTRA]bool isHidden() const)
	TEST_EQUAL(state_ptr->isHidden(), true)
	TEST_EQUAL(state_ptr2->isHidden(), false)
RESULT

CHECK([EXTRA]void setIdentifier(UInt id))
	state_ptr->setIdentifier(1234);
	TEST_EQUAL(state_ptr->getIdentifier(), 1234)
	state_ptr->setIdentifier(27);
RESULT

CHECK([EXTRA]void setHidden(bool hidden))
	state_ptr->setHidden(false);
	TEST_EQUAL(state_ptr->isHidden(), false)
	state_ptr->setHidden(true);
	TEST_EQUAL(state_ptr->isHidden(), true)
RESULT

CHECK([EXTRA]const std::set<HMMStateLight*>& getPredecessorStates() const)
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 0)
RESULT

CHECK([EXTRA]const std::set<HMMStateLight*>& getSuccessorStates() const)
	TEST_EQUAL(state_ptr->getSuccessorStates().size(), 0)
RESULT

CHECK([EXTRA]void addPredecessorState(HMMStateLight* state))
	state_ptr->addPredecessorState(state_ptr2);
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 1);
	TEST_EQUAL(*state_ptr->getPredecessorStates().begin(), state_ptr2);
RESULT

CHECK([EXTRA]void deletePredecessorState(HMMStateLight* state))
	state_ptr->deletePredecessorState(state_ptr2);
	TEST_EQUAL(state_ptr->getPredecessorStates().size(), 0);
RESULT

CHECK([EXTRA]void addSuccessorState(HMMStateLight* state))
  state_ptr->addSuccessorState(state_ptr2);
  TEST_EQUAL(state_ptr->getSuccessorStates().size(), 1);
  TEST_EQUAL(*state_ptr->getSuccessorStates().begin(), state_ptr2);
RESULT

CHECK([EXTRA]void deleteSuccessorState(HMMStateLight* state))
  state_ptr->deleteSuccessorState(state_ptr2);
  TEST_EQUAL(state_ptr->getSuccessorStates().size(), 0);
RESULT


// Hidden Markov Model Tests
CHECK(HiddenMarkovModelLight())
  ptr = new HiddenMarkovModelLight();
  TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(HiddenMarkovModelLight(const HiddenMarkovModelLight& hmm_new))
	HiddenMarkovModelLight copy(*ptr);
	TEST_EQUAL(copy.getNumberOfStates(), ptr->getNumberOfStates())
RESULT

CHECK(HiddenMarkovModelLight& operator = (const HiddenMarkovModelLight&))
  HiddenMarkovModelLight copy;
	copy = *ptr;
	TEST_EQUAL(copy.getNumberOfStates(), ptr->getNumberOfStates())
RESULT

CHECK(~HiddenMarkovModelLight())
  delete ptr;
RESULT

ptr = new HiddenMarkovModelLight();

CHECK(UInt getNumberOfStates() const)
  TEST_EQUAL(ptr->getNumberOfStates(), 0)
RESULT

CHECK(void addNewState(HMMStateLight* state))
  ptr->addNewState(state_ptr);
  TEST_EQUAL(ptr->getNumberOfStates(), 1)
  ptr->addNewState(state_ptr2);
  ptr->addNewState(state_ptr3);
RESULT

CHECK(HMMStateLight* getState(UInt id1))
  TEST_EQUAL(ptr->getState(27), state_ptr)
RESULT

CHECK(const HMMStateLight* getState(UInt id1) const)
	TEST_EQUAL(ptr->getState(27), state_ptr)
RESULT

CHECK(double getTransitionProbability(HMMStateLight*, HMMStateLight*) const)
  TEST_REAL_EQUAL(ptr->getTransitionProbability(27, 123), 0.0)
RESULT

CHECK(double getTransitionProbability(UInt id1, UInt id2) const)
  TEST_REAL_EQUAL(ptr->getTransitionProbability(state_ptr, state_ptr2), 0.0)
RESULT

CHECK(void setTransitionProbability(UInt id1, UInt id2, double prob))
  ptr->setTransitionProbability(27, 123, 0.3);
  TEST_REAL_EQUAL(ptr->getTransitionProbability(27, 123), 0.3)
RESULT

CHECK(void setTransitionProbability(HMMStateLight* s1, HMMStateLight* s2, double prob))
  ptr->setTransitionProbability(state_ptr, state_ptr2, 0.4);
  TEST_REAL_EQUAL(ptr->getTransitionProbability(state_ptr, state_ptr2), 0.4)
RESULT

CHECK(void addSynonymTransition(UInt name1, UInt name2, UInt synonym1, UInt synonym2))
  HMMStateLight* s1 = new HMMStateLight(28);
  HMMStateLight* s2 = new HMMStateLight(124);
  ptr->addNewState(s1);
  ptr->addNewState(s2);
  ptr->addSynonymTransition(27, 123, 28, 124);
RESULT

CHECK(void buildSynonyms())
  ptr->buildSynonyms();
  TEST_REAL_EQUAL(ptr->getTransitionProbability(28, 124), 0.4)
RESULT

CHECK(void setInitialTransitionProbability(UInt id, double prob))
RESULT

CHECK(void setInitialTransitionProbability(HMMStateLight* state, double prob))
RESULT

CHECK(void setTrainingEmissionProbability(UInt id, double prob))
RESULT

CHECK(void setTrainingEmissionProbability(HMMStateLight* state, double prob))
RESULT

CHECK(void enableTransition(HMMStateLight* s1, HMMStateLight* s2))
RESULT

CHECK(void enableTransition(UInt id1, UInt id2))
RESULT

CHECK(void disableTransition(HMMStateLight* s1, HMMStateLight* s2))
RESULT

CHECK(void disableTransition(UInt id1, UInt id2))
RESULT

CHECK(void disableTransitions())
  ptr->disableTransitions();
RESULT

CHECK(void calculateEmissionProbabilities(Map<HMMStateLight*, double>& emission_probs))
RESULT

CHECK(void train())
RESULT

CHECK(void evaluate())
RESULT

CHECK(void estimateUntrainedTransitions())
RESULT

CHECK(void clearInitialTransitionProbabilities())
RESULT

CHECK(void clearTrainingEmissionProbabilities())
RESULT

CHECK(void addIdToName(UInt id, const String& name))
RESULT

CHECK(void dump())
	// nothing to test
RESULT

CHECK(void forwardDump())
	// nothing to test
RESULT

CHECK(void write(std::ostream& out))
RESULT

CHECK(void writeGraphMLFile(const String& filename))
RESULT

CHECK(void readFromFile(const String& filename))
RESULT

CHECK(void setPseudoCounts(double pseudo_counts))
  ptr->setPseudoCounts(10e-3);
RESULT

CHECK(double getPseudoCounts() const)
  TEST_EQUAL(ptr->getPseudoCounts(), 10e-3)
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
