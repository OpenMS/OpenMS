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

// the HMM
// 0.5     0.25     0.25
//  |	      |        |
//  v       v        v
//  A       B        C
//  |       |        |
//  v       v        v
// A_2     B_2      C_2
//  \       /        |
//   \     /         |
//    \   /          |
//     v v           v
//    AB_3          C_3
//
// each edge is accompanied by an edge to the "end" state
// the weight of each edge is 0.5

HiddenMarkovModel the_hmm;
HMMState* state_A = new HMMState("A", true);
HMMState* state_B = new HMMState("B", true);
HMMState* state_C = new HMMState("C", true);
HMMState* state_A_2 = new HMMState("A_2", true);
HMMState* state_B_2 = new HMMState("B_2", true);
HMMState* state_C_2 = new HMMState("C_2", true);
HMMState* state_AB_3 = new HMMState("AB_3", false);
HMMState* state_C_3 = new HMMState("C_3", false);
HMMState* state_end = new HMMState("end", false);

the_hmm.addNewState(state_A);
the_hmm.addNewState(state_B);
the_hmm.addNewState(state_C);
the_hmm.addNewState(state_A_2);
the_hmm.addNewState(state_B_2);
the_hmm.addNewState(state_C_2);
the_hmm.addNewState(state_AB_3);
the_hmm.addNewState(state_C_3);
the_hmm.addNewState(state_end);

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
	TEST_EQUAL(ptr->getNumberOfStates(), 3)
END_SECTION

START_SECTION((HMMState* getState(const String& name)))
	TEST_EQUAL(ptr->getState("state_name_hidden"), state_ptr)
END_SECTION

START_SECTION((const HMMState* getState(const String& name) const))
	TEST_EQUAL(ptr->getState("state_name_hidden"), state_ptr)
END_SECTION

START_SECTION((DoubleReal getTransitionProbability(const String& s1, const String& s2) const))
	TEST_REAL_SIMILAR(ptr->getTransitionProbability("state_name_hidden", "state_name_emitting"), 0.0)
END_SECTION

START_SECTION((void setTransitionProbability(const String& s1, const String& s2, DoubleReal prob)))
	ptr->setTransitionProbability("state_name_hidden", "state_name_emitting", 0.3);
	TEST_REAL_SIMILAR(ptr->getTransitionProbability("state_name_hidden", "state_name_emitting"), 0.3)

	the_hmm.setTransitionProbability("A", "A_2", 0.5);
	the_hmm.setTransitionProbability("A", "end", 0.5);
	the_hmm.setTransitionProbability("B", "B_2", 0.5);
	the_hmm.setTransitionProbability("B", "end", 0.5);
	the_hmm.setTransitionProbability("C", "C_2", 0.5);
	the_hmm.setTransitionProbability("C", "end", 0.5);
	the_hmm.setTransitionProbability("A_2", "AB_3", 0.5);
	the_hmm.setTransitionProbability("A_2", "end", 0.5);
	the_hmm.setTransitionProbability("B_2", "AB_3", 0.5);
	the_hmm.setTransitionProbability("B_2", "end", 0.5);
	the_hmm.setTransitionProbability("C_2", "C_3", 0.5);
	the_hmm.setTransitionProbability("C_2", "end", 0.5);
	
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("A", "A_2"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("A", "end"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("B", "B_2"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("B", "end"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("C", "C_2"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("C", "end"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("A_2", "AB_3"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("A_2", "end"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("B_2", "AB_3"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("B_2", "end"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("C_2", "C_3"), 0.5)
	TEST_REAL_SIMILAR(the_hmm.getTransitionProbability("C_2", "end"), 0.5)

END_SECTION

START_SECTION((void addSynonymTransition(const String& name1, const String& name2, const String& synonym1, const String& synonym2)))
	HMMState* s1 = new HMMState("state_name_hidden2");
	HMMState* s2 = new HMMState("state_name_emitting2");
	ptr->addNewState(s1);
	ptr->addNewState(s2);
	ptr->addSynonymTransition("state_name_hidden", "state_name_emitting", "state_name_hidden2", "state_name_emitting2");
	NOT_TESTABLE
END_SECTION

START_SECTION((void setInitialTransitionProbability(const String& state, DoubleReal prob)))
	ptr->setInitialTransitionProbability("state_name_hidden2", 1.0);
	NOT_TESTABLE

	the_hmm.setInitialTransitionProbability("A", 0.5);
	the_hmm.setInitialTransitionProbability("B", 0.25);
	the_hmm.setInitialTransitionProbability("C", 0.25);

END_SECTION

START_SECTION((void enableTransition(const String& s1, const String& s2)))
	the_hmm.enableTransition("A", "A_2");
	the_hmm.enableTransition("A", "end");
	the_hmm.enableTransition("B", "B_2");
	the_hmm.enableTransition("B", "end");
	the_hmm.enableTransition("C", "C_2");
	the_hmm.enableTransition("C", "end");
	the_hmm.enableTransition("A_2", "AB_3");
	the_hmm.enableTransition("A_2", "end");
	the_hmm.enableTransition("B_2", "AB_3");
	the_hmm.enableTransition("B_2", "end");
	the_hmm.enableTransition("C_2", "C_3");
	the_hmm.enableTransition("C_2", "end");
	NOT_TESTABLE // will be tested implicetly below
END_SECTION

START_SECTION((void disableTransition(const String& s1, const String& s2)))
	NOT_TESTABLE // will be tested implicitely below
END_SECTION

START_SECTION((void disableTransitions()))
	ptr->disableTransitions();
	NOT_TESTABLE
END_SECTION

START_SECTION((void calculateEmissionProbabilities(Map<HMMState*, DoubleReal>& emission_probs)))
	Map<HMMState*, DoubleReal> emission_probs;
	the_hmm.calculateEmissionProbabilities(emission_probs);
	TEST_EQUAL(emission_probs.size(), 3)
	DoubleReal sum(0);
	TOLERANCE_ABSOLUTE(0.01)
	for (Map<HMMState*, DoubleReal>::ConstIterator it = emission_probs.begin(); it != emission_probs.end(); ++it)
	{
		if (it->first->getName() == "end")
		{
			sum += it->second;
			TEST_REAL_SIMILAR(it->second, 12.0/16.0)
		}
		else if (it->first->getName() == "AB_3")
		{
			sum += it->second;
			TEST_REAL_SIMILAR(it->second, 3.0/16.0)
		}
		else if (it->first->getName() == "C_3")
		{
			sum += it->second;
			TEST_REAL_SIMILAR(it->second, 1.0/16.0)
		}
	}
	TEST_REAL_SIMILAR(sum , 1.0)
END_SECTION

START_SECTION((void setTrainingEmissionProbability(const String& state, DoubleReal prob)))
	the_hmm.setTrainingEmissionProbability("end", 0.5);
	the_hmm.setTrainingEmissionProbability("AB_3", 0.3);
	the_hmm.setTrainingEmissionProbability("C_3", 0.2);
	NOT_TESTABLE
END_SECTION

START_SECTION((void train()))
	the_hmm.train();
	NOT_TESTABLE
END_SECTION

START_SECTION((void evaluate()))
	the_hmm.evaluate();
	NOT_TESTABLE
END_SECTION

START_SECTION((void estimateUntrainedTransitions()))
	NOT_TESTABLE // only applicable to the fragmentation model
END_SECTION

START_SECTION(([EXTRA] void calculateEmissionProbabilities(Map<HMMState*, DoubleReal>& emission_probs)))
	Map<HMMState*, DoubleReal> emission_probs;
  the_hmm.calculateEmissionProbabilities(emission_probs);
  TEST_EQUAL(emission_probs.size(), 3)
  DoubleReal sum(0);
	TOLERANCE_ABSOLUTE(0.01)
  for (Map<HMMState*, DoubleReal>::ConstIterator it = emission_probs.begin(); it != emission_probs.end(); ++it)
  {
    if (it->first->getName() == "end")
    {
      sum += it->second;
      TEST_REAL_SIMILAR(it->second, 0.8456)
    }
    else if (it->first->getName() == "AB_3")
    {
      sum += it->second;
      TEST_REAL_SIMILAR(it->second, 0.125)
    }
    else if (it->first->getName() == "C_3")
    {
      sum += it->second;
      TEST_REAL_SIMILAR(it->second, 0.02941)
    }
  }
  TEST_REAL_SIMILAR(sum , 1.0)

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
	NOT_TESTABLE
END_SECTION

START_SECTION((void clearTrainingEmissionProbabilities()))
	ptr->clearTrainingEmissionProbabilities();
	NOT_TESTABLE
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
	String str_stream = ss.str();
	TEST_EQUAL(str_stream.hasSubstring("State"), true);
	TEST_EQUAL(str_stream.hasSubstring("Transition"), true);
	TEST_EQUAL(str_stream.hasSubstring("Synonym"), true);
END_SECTION

START_SECTION((void writeGraphMLFile(const String& filename)))
	String filename;
	NEW_TMP_FILE(filename)
	ptr->writeGraphMLFile(filename);
	//TEST_FILE_SIMILAR(filename, OPENMS_GET_TEST_DATA_PATH("HiddenMarkovModel_test.graphML"))
	NOT_TESTABLE // just a convenience function; the sorting of the nodes will depend on the instance...
END_SECTION

START_SECTION((void setVariableModifications(const StringList &modifications)))
	StringList mods = StringList::create("Carboxymethyl (C),Oxidation (M)");
	ptr->setVariableModifications(mods);
	NOT_TESTABLE
END_SECTION

START_SECTION((void clear()))
	ptr->clear();
	TEST_EQUAL(ptr->getNumberOfStates(), 0)
END_SECTION

START_SECTION(void addNewState(const String &name))
	ptr->addNewState("new_fancy_state");
	TEST_EQUAL(ptr->getNumberOfStates(), 1)
END_SECTION
  
START_SECTION(void setPseudoCounts(DoubleReal pseudo_counts))
	ptr->setPseudoCounts(10e-3);
	NOT_TESTABLE // tested in next section
END_SECTION

START_SECTION(DoubleReal getPseudoCounts() const)
	TEST_EQUAL(ptr->getPseudoCounts(), 10e-3)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST

