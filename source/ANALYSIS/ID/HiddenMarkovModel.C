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


#include <OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>

#include <iostream>
#include <fstream>
#include <stack>

//#define PSEUDO_COUNTS 0.0000001
#define PSEUDO_COUNTS   1e-15

#define HIDDEN_MARKOV_MODEL_DEBUG
#undef HIDDEN_MARKOV_MODEL_DEBUG

using namespace std;

namespace OpenMS 
{
	// HMMState implementation
	HMMState::HMMState()
		: hidden_(true)
	{
	}

	HMMState::HMMState(const String& name, bool hidden)
		: hidden_(hidden),
			name_(name)
	{
	}

	HMMState::HMMState(const HMMState& state)
		: hidden_(state.hidden_),
			name_(state.name_),
			pre_states_(state.pre_states_),
			succ_states_(state.succ_states_)
	{
	}
	
	HMMState::~HMMState()
	{
	}
	
	void HMMState::setName(const String& name)
	{
		name_ = name;
	}

	const String& HMMState::getName() const
	{
		return name_;
	}

	HMMState& HMMState::operator = (const HMMState& state)
	{
		hidden_ = state.hidden_;
		name_ = state.name_;
		pre_states_ = state.pre_states_;
		succ_states_ = state.succ_states_;
		return *this;
	}

	void HMMState::addSuccessorState(HMMState* state)
	{
		succ_states_.insert(state);
	}

	void HMMState::deleteSuccessorState(HMMState* state)
	{
		succ_states_.erase(state);
	}

	void HMMState::addPredecessorState(HMMState* state)
	{
		pre_states_.insert(state);
	}

	void HMMState::deletePredecessorState(HMMState* state)
	{
		pre_states_.erase(state);
	}

	const set<HMMState::HMMState*>& HMMState::getPredecessorStates() const
	{
		return pre_states_;
	}

	const set<HMMState::HMMState*>& HMMState::getSuccessorStates() const
	{
		return succ_states_;
	}

	bool HMMState::isHidden() const
	{
		return hidden_;
	}

	void HMMState::setHidden(bool hidden)
	{
		hidden_ = hidden;
	}

	// The hidden markov model
	HiddenMarkovModel::HiddenMarkovModel() 
	{
	}

	HiddenMarkovModel::HiddenMarkovModel(const HiddenMarkovModel& hmm)
		:	trans_(hmm.trans_),
			count_trans_(hmm.count_trans_),
			name_to_state_(hmm.name_to_state_),
			states_(hmm.states_)
	{
	}
	
	HiddenMarkovModel::~HiddenMarkovModel()
	{
	}

	HiddenMarkovModel& HiddenMarkovModel::operator = (const HiddenMarkovModel& hmm)
	{
		trans_ = hmm.trans_;
		count_trans_ = hmm.count_trans_;
		name_to_state_ = hmm.name_to_state_;
		states_ = hmm.states_;
		return *this;
	}

	HMMState* HiddenMarkovModel::getState(const String& name)
	{
		if (!name_to_state_.has(name))
		{
			throw Exception::ElementNotFound<String>(__FILE__, __LINE__, __PRETTY_FUNCTION__, name);
		}
		return name_to_state_[name];
	}

	const HMMState* HiddenMarkovModel::getState(const String& name) const
	{
		if (!name_to_state_.has(name))
		{
			throw Exception::ElementNotFound<String>(__FILE__, __LINE__, __PRETTY_FUNCTION__, name);
		}
		return name_to_state_[name];
	}
	
	void HiddenMarkovModel::writetoYGFFile(const String& filename)
	{
		set<HMMState*> states;
		HashMap<HMMState*, vector<HMMState*> > transitions;
		
		stack<HMMState*> s;
		for (HashMap<HMMState*, double>::ConstIterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
		{
			s.push(it->first);
		}

		while (s.size() != 0)
		{
			HMMState* state = s.top();
			states.insert(state);
			s.pop();
			for (set<HMMState*>::const_iterator it = state->getSuccessorStates().begin(); it != state->getSuccessorStates().end(); ++it)
			{
				transitions[state].push_back(*it);
				s.push(*it);
			}
		}

	
		ofstream out(filename.c_str());
    out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << endl;

		out << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns/graphml\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" " <<
         "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns/graphml http://www.yworks.com/xml/schema/graphml/1.0/ygraphml.xsd\" " <<
				 "xmlns:y=\"http://www.yworks.com/xml/graphml\">" << endl;

		out << "<key id=\"d0\" for=\"node\" yfiles.type=\"nodegraphics\"/>" << endl;
  	out << "<key id=\"d1\" for=\"edge\" yfiles.type=\"edgegraphics\"/>" << endl;
    out << "  <graph id=\"G\" edgedefault=\"directed\">" << endl;
		for (set<HMMState*>::const_iterator it = states.begin(); it != states.end(); ++it)
		{
			out << "    <node id=\"" << (*it)->getName() << "\">" << endl;
			out << "      <data key=\"d0\">" << endl;
			out << "        <y:ShapeNode>" << endl;
			out << "          <y:NodeLabel>" << (*it)->getName() << "</y:NodeLabel>" << endl;
			out << "        </y:ShapeNode>" << endl;
			out << "      </data>" << endl;
			out << "    </node>" << endl;
		}

		for (HashMap<HMMState*, vector<HMMState*> >::ConstIterator it = transitions.begin(); it != transitions.end(); ++it)
		{
			for (vector<HMMState*>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				out << "    <edge source=\"" << it->first->getName() << "\" target=\"" << (*it1)->getName() << "\">" << endl;
				out << "      <data key=\"d1\">" << endl;
				out << "        <y:PolyLineEdge>" << endl;
				//out << "          <y:EdgeLabel>" << train_count_trans_[it->first][*it1] << "</y:EdgeLabel>" << endl;
				//cerr << it->first->getName() << " " << (*it1)->getName() << " " << train_count_trans_[it->first][*it1] << endl;
				out << "          <y:EdgeLabel>" << getTransitionProbability(it->first, *it1) << "</y:EdgeLabel>" << endl;
				out << "        </y:PolyLineEdge>" << endl;
				out << "      </data>" << endl;
				out << "    </edge>" << endl;
			}
		}

		out << "  </graph>" << endl;
		out << "</graphml>" << endl;
		
	}

	void HiddenMarkovModel::write(ostream& out)
	{
		// write states
		for (set<HMMState*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
		{
			out << "State " << (*it)->getName();

			if (!(*it)->isHidden())
			{
				out << " false";
			}
			out << endl;
		}

		// write transitions
		for (HashMap<HMMState*, HashMap<HMMState*, double> >::ConstIterator it1 = trans_.begin(); it1 != trans_.end(); ++it1)
		{
			for (HashMap<HMMState*, double>::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				out << "Transition " << it1->first->getName() << " " << it2->first->getName() << " " << it2->second << endl;
			}
		}
	
		// write synonyms
		/*
		for (HashMap<String, HashMap<String, pair<String, String> > >::ConstIterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
		{
			for (HashMap<String, pair<String, String> >::ConstIterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				out << "Synonym " << it->first << " " << it2->first << " " << it2->second.first << " " << it2->second.second << endl;
			}
		}*/

    for (HashMap<HMMState*, HashMap<HMMState*, std::pair<HMMState*, HMMState*> > >::ConstIterator it1 = synonym_trans_.begin(); it1 != synonym_trans_.end(); ++it1)
    {
      for (HashMap<HMMState*, std::pair<HMMState*, HMMState*> >::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
				//cerr << "Synonym " << it1->first->getName() << " " << it2->first->getName() << " " << it2->second.first->getName() << " " << it2->second.second->getName() << endl;
        out << "Synonym " << it1->first->getName() << " " << it2->first->getName() << " " << it2->second.first->getName() << " " << it2->second.second->getName() << endl;
      }
    }

		return;
	}

	void HiddenMarkovModel::clear()
	{
		for (set<HMMState*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
    {
      delete *it;
    }
    trans_.clear();
    count_trans_.clear();
    forward_.clear();
    backward_.clear();
    name_to_state_.clear();
    train_emission_prob_.clear();
    init_prob_.clear();
    states_.clear();
    trained_trans_.clear();
    synonym_trans_.clear();
    synonym_trans_names_.clear();	
		return;
	}
	
	/*
	void HiddenMarkovModel::readFromFile(const String& filename)
	{
		// clear old model
		for (set<HMMState*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
		{
			delete *it;
		}
		trans_.clear();
		count_trans_.clear();
		forward_.clear();
		backward_.clear();
		name_to_state_.clear();
		train_emission_prob_.clear();
		init_prob_.clear();
		states_.clear();
		trained_trans_.clear();
		synonym_trans_.clear();
		synonym_trans_names_.clear();
		
		ifstream in(filename.c_str());
		char buf[10000];
		HMMState* s = 0;
		HashMap<HMMState*, vector<String> > pre_states;
		HashMap<HMMState*, vector<String> > succ_states;
		while (in.getline(buf, 10000, '\n'))
		{
			String line(buf);
			if (line == "#trans_:")
			{
				break;
			}
			if (line == "#states_:")
			{
				continue;
			}
			vector<String> split;
			line.split(' ', split);
			String split_0;
			if (split.size() == 0)
			{
				split_0 = line;
			}
			else
			{
				split_0 = split[0];
			}
			if (split_0 == String("name_:"))
			{
				s = new HMMState(split[1]);
			}
			if (split_0 == "hidden_:")
			{
				if (split[1] == "false")
				{
					s->setHidden(false);
				}
			
				states_.insert(s);
				name_to_state_[s->getName()] = s;
			}
		}

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "states block ready (#states=" << states_.size() << ", name_to_state_.size()=" << name_to_state_.size() << ")" << endl;
		#endif

		// trans_ block
		while (in.getline(buf, 10000, '\n'))
		{
			String line(buf);
			line.trim();
			if (line[0] == '#' && line == "#synonym_trans_:")
			{
				break;
			}
			vector<String> split;
			line.split(' ', split);
			trans_[name_to_state_[split[0]]][name_to_state_[split[2]]] = atof(split[3].c_str());
		}


		// synonym trans block
		while (in.getline(buf, 10000, '\n'))
		{
			String line(buf);
			line.trim();
			vector<String> split;
			line.split(' ', split);
			synonym_trans_names_[split[0]][split[1]] = make_pair<String, String>(split[2], split[3]);
		}

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "trans_.size()=" << trans_.size() << ", synonym_trans_names_.size()=" << synonym_trans_names_.size() << endl;
		#endif

		buildSynonyms();

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		dump();
		#endif
	}
*/
	
	Size HiddenMarkovModel::getNumberOfStates() const
	{
		return states_.size();
	}

	void HiddenMarkovModel::addNewState(HMMState* s)
	{
		states_.insert(s);
		if (!name_to_state_.has(s->getName()))
		{
			name_to_state_[s->getName()] = s;
		}
		else
		{
			cerr << "HiddenMarkovModel: state name '" << s->getName() << "' (" << s << ") already used!" << endl;
		}
	}

	void HiddenMarkovModel::addSynonymTransition(const String& name1, const String& name2, const String& synonym1, const String& synonym2)
	{
		if (!name_to_state_.has(name1))
		{
			cerr << "state '" << name1 << "' unknown" << endl;
		}
    if (!name_to_state_.has(name2))
    {
      cerr << "state '" << name2 << "' unknown" << endl;
    }
    if (!name_to_state_.has(synonym1))
    {
      cerr << "state '" << synonym1 << "' unknown" << endl;
    }
    if (!name_to_state_.has(synonym2))
    {
      cerr << "state '" << synonym2 << "' unknown" << endl;
    }
		synonym_trans_names_[synonym1][synonym2] = make_pair<String, String>(name1, name2);
	}

	void HiddenMarkovModel::setTransitionProbability(HMMState * s1, HMMState * s2, double trans_prob)
	{
		trans_[s1][s2] = trans_prob;
		s1->addSuccessorState(s2);
		s2->addPredecessorState(s1);
		enabled_trans_[s1].insert(s2);
		training_steps_count_[s1][s2] = 0;
	}

	void HiddenMarkovModel::setTransitionProbability(const String& s1, const String& s2, double trans_prob)
	{
		trans_[name_to_state_[s1]][name_to_state_[s2]] = trans_prob;
		name_to_state_[s1]->addSuccessorState(name_to_state_[s2]);
		name_to_state_[s2]->addPredecessorState(name_to_state_[s1]);
		enabled_trans_[name_to_state_[s1]].insert(name_to_state_[s2]);
		training_steps_count_[name_to_state_[s1]][name_to_state_[s2]] = 0;
	}

	double HiddenMarkovModel::getTransitionProbability(const String& s1, const String& s2) const
	{
		if (!name_to_state_.has(s1))
		{
			throw Exception::ElementNotFound<String>(__FILE__, __LINE__, __PRETTY_FUNCTION__, s1);
		}
		if (!name_to_state_.has(s2))
		{
			throw Exception::ElementNotFound<String>(__FILE__, __LINE__, __PRETTY_FUNCTION__, s2);
		}
		return getTransitionProbability(name_to_state_[s1], name_to_state_[s2]);
	}
	
	double HiddenMarkovModel::getTransitionProbability(HMMState* s1, HMMState* s2) const
	{
		HMMState* state1 = s1;
		HMMState* state2 = s2;
		
		// check if transition is a synonym transition
		if (synonym_trans_.has(state1) && synonym_trans_[state1].has(state2))
		{
			HMMState* tmp = synonym_trans_[state1][state2].first;
			state2 = synonym_trans_[state1][state2].second; //name_to_state_[p.second];
			state1 = tmp;
		}	
		
		// get the transition prob 
		if (trans_.has(state1) && trans_[state1].has(state2))
		{
			return trans_[state1][state2];
		}
		else
		{
			return 0;
		}
	}

	void HiddenMarkovModel::train()
	{
		trained_trans_.clear();
		train_count_trans_.clear();

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "calc forward part" << endl;
		#endif
		
		calculateForwardPart_();

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "calc backward part" << endl;
		#endif
		
		calculateBackwardPart_();
	

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "calc px" << endl;
		#endif
	
		// calc p_x from forward part
		double px(0);

		for (HashMap<HMMState*, double>::ConstIterator it1 = train_emission_prob_.begin(); it1 != train_emission_prob_.end(); ++it1)
		{
			for (set<HMMState*>::const_iterator it2 = it1->first->getPredecessorStates().begin(); it2 != it1->first->getPredecessorStates().end(); ++it2)
			{
				px += getForwardVariable_(*it2);
			}
		}

		double num_px(0);
		if (px != 0)
		{
			num_px = 1/px;
		}

#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "px=" << num_px << endl;
		cerr << "add contributions to count_trans" << endl;
#endif

		// add contributions to count_trans_
		for (set<pair<HMMState*, HMMState*> >::const_iterator it = trained_trans_.begin(); it != trained_trans_.end(); ++it)
		{
			double tmp(0);
			tmp = num_px * getForwardVariable_(it->first) * getBackwardVariable_(it->second) * getTransitionProbability(it->first, it->second);
			tmp += PSEUDO_COUNTS;
			HMMState* s1 = it->first;
			HMMState* s2 = it->second;
			train_count_trans_[s1][s2] = tmp;

			if (synonym_trans_.has(s1) && synonym_trans_[s1].has(s2))
			{
				HMMState* tmp = synonym_trans_[s1][s2].first; //name_to_state_[p.first];
				s2 = synonym_trans_[s1][s2].second; // name_to_state_[p.second];
				s1 = tmp;
			}
			
			if (count_trans_.has(s1) && count_trans_[s1].has(s2))
			{
				count_trans_[s1][s2] += tmp;
			}
			else
			{
				count_trans_[s1][s2] = tmp;
			}
			training_steps_count_[s1][s2]++;
		}
	}

	void HiddenMarkovModel::calculateForwardPart_()
	{
		forward_.clear();
		set<HMMState*> succ;
		for (HashMap<HMMState*, double>::Iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
		{
			forward_[it->first] = it->second;
		}
		
		for (HashMap<HMMState*, double>::Iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
		{
			succ.insert(it->first->getSuccessorStates().begin(), it->first->getSuccessorStates().end());
		
			while (succ.size() != 0)
			{
				set<HMMState*> succ_new;
				for (set<HMMState*>::const_iterator it = succ.begin(); it != succ.end(); ++it)
				{
					set<HMMState*> pre = (*it)->getPredecessorStates();
					double sum(0);
					for (set<HMMState*>::const_iterator it2 = pre.begin(); it2 != pre.end(); ++it2)
					{
						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
						cerr << "\tadding from " << (*it2)->getName() << " " << getForwardVariable_(*it2) << " * " << getTransitionProbability(*it2, *it) << endl;
						#endif
						sum += getForwardVariable_(*it2) * getTransitionProbability(*it2, *it);
						trained_trans_.insert(make_pair<HMMState*, HMMState*>(*it2, *it));
					}
					forward_[*it] = sum;
					succ_new.insert((*it)->getSuccessorStates().begin(), (*it)->getSuccessorStates().end());
					#ifdef HIDDEN_MARKOV_MODEL_DEBUG
					cerr << "f: " << (*it)->getName() << "\t" << sum << endl;
					#endif
				}
				succ = succ_new;
			}
		}
	}

	void HiddenMarkovModel::calculateBackwardPart_()
	{
		backward_.clear();
		set<HMMState*> pre;
		for (HashMap<HMMState*, double>::Iterator it = train_emission_prob_.begin(); it != train_emission_prob_.end(); ++it)
		{
			backward_[it->first] = it->second;
		}

		for (HashMap<HMMState*, double>::Iterator it = train_emission_prob_.begin(); it != train_emission_prob_.end(); ++it)
		{
			#ifdef HIDDEN_MARKOV_MODEL_DEBUG
			cerr << "b:" << it->first->getName() << " " << it->second << endl;
			#endif
			pre.insert(it->first->getPredecessorStates().begin(), it->first->getPredecessorStates().end());
		
			while (pre.size() != 0)
			{
				set <HMMState*> pre_new;
				for (set<HMMState*>::const_iterator it = pre.begin(); it != pre.end(); ++it)
				{
					set<HMMState*> succ = (*it)->getSuccessorStates();
					double sum(0);
					for (set<HMMState*>::const_iterator it2 = succ.begin(); it2 != succ.end(); ++it2)
					{
						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
						cerr << "\tadding from " << (*it2)->getName() << " " << getBackwardVariable_(*it2) << " * " << getTransitionProbability(*it, *it2) << endl;
						#endif
						sum += getBackwardVariable_(*it2) * getTransitionProbability(*it, *it2);
						trained_trans_.insert(make_pair<HMMState*, HMMState*>(*it, *it2));
					}
					backward_[*it] = sum;
					pre_new.insert((*it)->getPredecessorStates().begin(), (*it)->getPredecessorStates().end());
					#ifdef HIDDEN_MARKOV_MODEL_DEBUG
					cerr << "b: " << (*it)->getName() << "\t" << sum << endl;
					#endif
				}
				pre = pre_new;
			}
		}
	}

	double HiddenMarkovModel::getForwardVariable_(HMMState* state)
	{
		return forward_.has(state) ? forward_[state] : 0;
	}

	double HiddenMarkovModel::getBackwardVariable_(HMMState* state)
	{
		return backward_.has(state) ? backward_[state] : 0;
	}

	void HiddenMarkovModel::evaluate()
	{
		for (HashMap<HMMState*, HashMap<HMMState*, double> >::ConstIterator it1 = count_trans_.begin(); it1 != count_trans_.end(); ++it1)
		{
			double sum(0);
			for (HashMap<HMMState*, double>::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				if (count_trans_.has(it1->first) && count_trans_[it1->first].has(it2->first))
				{
					sum += count_trans_[it1->first][it2->first];
				}
			}

			if (sum != 0)
			{
				for (HashMap<HMMState*, double>::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
				{
					if (count_trans_.has(it1->first) && count_trans_[it1->first].has(it2->first))
					{
						trans_[it1->first][it2->first] = count_trans_[it1->first][it2->first] / sum;
					}
				}
			}
		}
	}

	void HiddenMarkovModel::setInitialTransitionProbability(const String& state, double prob)
	{
		init_prob_[name_to_state_[state]] = prob;
	}

	void HiddenMarkovModel::clearInitialTransitionProbabilities()
	{
		init_prob_.clear();
	}

	void HiddenMarkovModel::setTrainingEmissionProbability(const String& state, double prob)
	{
		train_emission_prob_[name_to_state_[state]] = prob;
	}

	void HiddenMarkovModel::clearTrainingEmissionProbabilities()
	{
		train_emission_prob_.clear();
	}

	void HiddenMarkovModel::enableTransition(const String& s1, const String& s2)
	{
		enableTransition(name_to_state_[s1], name_to_state_[s2]);
	}

	void HiddenMarkovModel::enableTransition(HMMState* s1, HMMState* s2)
	{
		s1->addSuccessorState(s2);
		s2->addPredecessorState(s1);
		enabled_trans_[s1].insert(s2);
	}

	void HiddenMarkovModel::disableTransition(const String& s1, const String& s2)
	{
		disableTransition(name_to_state_[s1], name_to_state_[s2]);
	}
	
	void HiddenMarkovModel::disableTransition(HMMState* s1, HMMState* s2)
	{
		s1->deleteSuccessorState(s2);
		s2->deletePredecessorState(s1);
		enabled_trans_[s1].erase(s2);
	}

	void HiddenMarkovModel::disableTransitions()
	{
		for (HashMap<HMMState*, set<HMMState*> >::ConstIterator it = enabled_trans_.begin(); it != enabled_trans_.end(); ++it)
		{
			for (set<HMMState*>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				disableTransition(it->first, *it2);
			}
		}
		enabled_trans_.clear();
	}

	void HiddenMarkovModel::calculateEmissionProbabilities(HashMap<HMMState*, double>& emission_probs)
	{
		HashMap<HMMState*, double> states = init_prob_;

		while (states.size() != 0)
		{
			HashMap<HMMState*, double> tmp = states;
			for (HashMap<HMMState*, double>::ConstIterator it = tmp.begin(); it != tmp.end(); ++it)
			{
				for (set<HMMState*>::const_iterator it2 = it->first->getSuccessorStates().begin(); it2 != it->first->getSuccessorStates().end(); ++it2)
				{
					if (states.has(*it2))
					{
						states[*it2] += states[it->first] * getTransitionProbability(it->first, *it2);
					}
					else
					{
						states[*it2] = states[it->first] * getTransitionProbability(it->first, *it2);
					}
					if (!(*it2)->isHidden())
					{
						if (emission_probs.has(*it2))
						{
							emission_probs[*it2] += states[*it2];
						}
						else
						{
							emission_probs[*it2] = states[*it2];
						}
					}
				}
				states.erase(it->first);
			}
		}
		
		return;
	}

	void HiddenMarkovModel::dump()
	{
		cerr << "dump of transitions: " << endl;
		for (HashMap<HMMState*, HashMap<HMMState*, double> >::ConstIterator it = trans_.begin(); it != trans_.end(); ++it)
		{
			for (HashMap<HMMState*, double>::ConstIterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				cout << it->first->getName() << " -> " << it1->first->getName() << " " << it1->second << endl;
			}
		}
		
		cerr << "dump completed" << endl;
	}

	void HiddenMarkovModel::forwardDump()
	{
   set<HMMState*> succ;
    for (HashMap<HMMState*, double>::Iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
    {
      succ.insert(it->first->getSuccessorStates().begin(), it->first->getSuccessorStates().end());

      while (succ.size() != 0)
      {
        set<HMMState*> succ_new;
        for (set<HMMState*>::const_iterator it = succ.begin(); it != succ.end(); ++it)
        {
					cerr << (*it)->getName() << endl;
          succ_new.insert((*it)->getSuccessorStates().begin(), (*it)->getSuccessorStates().end());
        }
        succ = succ_new;
      }
    }
	}

	void HiddenMarkovModel::buildSynonyms()
	{
		for (HashMap<String, HashMap<String, std::pair<String, String> > >::ConstIterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
		{
			for (HashMap<String, pair<String, String> >::ConstIterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				synonym_trans_[name_to_state_[it->first]][name_to_state_[it2->first]] = 
						make_pair<HMMState*, HMMState*>(name_to_state_[it2->second.first], name_to_state_[it2->second.second]);
			}
		}
	}

	void HiddenMarkovModel::estimateUntrainedTransitions()
	{
		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "Transition training stats" << endl;
		for (HashMap<HMMState*, HashMap<HMMState*, Size> >::ConstIterator it = training_steps_count_.begin(); it != training_steps_count_.end(); ++it)
		{
			for (HashMap<HMMState*, Size>::ConstIterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				cerr << it->first->getName() << " " << it1->first->getName() << " " << it1->second << endl;
			}
		}

		cerr << "Estimation" << endl;
		#endif
		String residues("ACDEFGHIKMNPQRSTVWY");
		
		vector<String> suffixe;
		HMMState* s2 = name_to_state_["bxyz"];
		HMMState* end_state = name_to_state_["end"];
		for (Size i = 0; i != residues.size(); ++i)
		{
			s2 = name_to_state_["bxyz"];
			for (Size j = 0; j != residues.size(); ++j)
			{
				String aa1(residues[i]), aa2(residues[j]);
				if (training_steps_count_[name_to_state_[aa1 + aa2 + "bxyz"]][s2] == 0)
				{
					Size count(0);
					double sum(0);
					for (Size k = 0; k != residues.size(); ++k)
					{
						Size tmp = training_steps_count_[name_to_state_[aa1 + residues[k] + "bxyz"]][s2];
						if (tmp != 0)
						{
							sum += trans_[name_to_state_[aa1 + residues[k] + "bxyz"]][s2];
							count++;
						}
					}
					for (Size k = 0; k != residues.size(); ++k)
					{
						Size tmp = training_steps_count_[name_to_state_[residues[k] + aa2 +"bxyz"]][s2];
						if (tmp != 0)
						{
							sum += trans_[name_to_state_[residues[k] + aa2 +"bxyz"]][s2];
							count++;
						}
					}

					if (count != 0)
					{
						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
						cerr << "setting transitions of " << aa1 << aa2 << "bxyz -> bxyz to " << sum/double(count) << endl;
						#endif
						trans_[name_to_state_[aa1 + aa2 + "bxyz"]][s2] = sum/double(count);
						trans_[name_to_state_[aa1 + aa2 + "bxyz"]][end_state] = 1 - sum/double(count);
					}
				}
			}

			s2 = name_to_state_["axyz"];
      for (Size j = 0; j != residues.size(); ++j)
      {
        String aa1(residues[i]), aa2(residues[j]);
        if (training_steps_count_[name_to_state_[aa1 + aa2 + "axyz"]][s2] == 0)
        {
          Size count(0);
          double sum(0);
          for (Size k = 0; k != residues.size(); ++k)
          {
            Size tmp = training_steps_count_[name_to_state_[aa1 + residues[k] + "axyz"]][s2];
            if (tmp != 0)
            {
              sum += trans_[name_to_state_[aa1 + residues[k] + "axyz"]][s2];
              count++;
            }
          }
          for (Size k = 0; k != residues.size(); ++k)
          {
            Size tmp = training_steps_count_[name_to_state_[residues[k] + aa2 +"axyz"]][s2];
            if (tmp != 0)
            {
              sum += trans_[name_to_state_[residues[k] + aa2 +"axyz"]][s2];
              count++;
            }
          }
					if (count != 0)
					{
						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
	          cerr << "setting transitions of " << aa1 << aa2 << "axyz -> axyz to " << sum/double(count) << endl;
						#endif
  	        trans_[name_to_state_[aa1 + aa2 + "axyz"]][s2] = sum/double(count);
    	      trans_[name_to_state_[aa1 + aa2 + "axyz"]][end_state] = 1 - sum/double(count);
					}
        }
      }

			// sc and cr
			String sc_residues("HKDE");
		
			for (Size j = 0; j != sc_residues.size(); ++j)
			{
				String aa1(residues[i]), sc_res(sc_residues[j]);
				s2 = name_to_state_[sc_res];
				if (training_steps_count_[name_to_state_[aa1 + sc_res]][s2] == 0)
				{
					Size count(0);
					double sum(0);
					for (Size k = 0; k != residues.size(); ++k)
					{
						HMMState* s1 = name_to_state_[residues[k] + sc_res];
						Size tmp = training_steps_count_[s1][s2];
						if (tmp != 0)
						{
							sum += trans_[s1][s2];
							count++;
						}
					}

					if (count != 0)
					{
						trans_[name_to_state_[aa1 + sc_res]][s2] = sum/double(count);
						trans_[name_to_state_[aa1 + sc_res]][end_state] = 1 - sum/double(count);
					}
				}
			}

			String aa1(residues[i]), sc_res("RSC");
			s2 = name_to_state_["R"];

      if (training_steps_count_[name_to_state_[aa1 + sc_res]][s2] == 0)
      {
        Size count(0);
        double sum(0);
        for (Size k = 0; k != residues.size(); ++k)
        {
					HMMState* s1 = name_to_state_[residues[k] + sc_res];
          Size tmp = training_steps_count_[s1][s2];
          if (tmp != 0)
          {
            sum += trans_[s1][s2];
            count++;
          }
        }
				if (count != 0)
				{
	        trans_[name_to_state_[aa1 + sc_res]][s2] = sum/double(count);
	      	trans_[name_to_state_[aa1 + sc_res]][end_state] = 1 - sum/double(count);
				}
			}
		}

	}
}


