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

// TODO: !!!! niemals eine struktur aendern ueber die man iteriert!

#include <OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>

#include <iostream>
#include <fstream>
#include <stack>
#include <cmath>
#include <numeric>

#define SIMPLE_DEBUG2
#undef  SIMPLE_DEBUG2

#define HIDDEN_MARKOV_MODEL_DEBUG
#undef HIDDEN_MARKOV_MODEL_DEBUG

#define STATE_DEBUG
#undef STATE_DEBUG

#define EVALUATE_DEBUG
#undef EVALUATE_DEBUG

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
		#ifdef STATE_DEBUG
		cerr << "State: " << name << ", hidden=" << hidden << endl;
		#endif
	}

	HMMState::HMMState(const HMMState& state)
		: hidden_(state.hidden_),
			name_(state.name_)
	{
	}
	/*		pre_states_(state.pre_states_),
			succ_states_(state.succ_states_)
	{
		
	}*/
	
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
		pre_states_.clear();
		succ_states_.clear();
		//pre_states_ = state.pre_states_;
		//succ_states_ = state.succ_states_;
		return *this;
	}

	void HMMState::addSuccessorState(HMMState* state)
	{
#ifdef SIMPLE_DEBUG
		cerr << "'" << state->getName() << "'" << endl;
#endif
		succ_states_.insert(state);
	}

	void HMMState::deleteSuccessorState(HMMState* state)
	{
#ifdef SIMPLE_DEBUG
    cerr << "'" << state->getName() << "'" << endl;
#endif
		succ_states_.erase(state);
	}

	void HMMState::addPredecessorState(HMMState* state)
	{
#ifdef SIMPLE_DEBUG
    cerr << "'" << state->getName() << "'" << endl;
#endif
		pre_states_.insert(state);
	}

	void HMMState::deletePredecessorState(HMMState* state)
	{
#ifdef SIMPLE_DEBUG
    cerr << "'" << state->getName() << "'" << endl;
#endif
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
		: pseudo_counts_(0.0)
	{
	}

	HiddenMarkovModel::HiddenMarkovModel(const HiddenMarkovModel& hmm)
	{
		copy_(hmm);
	}
	
	HiddenMarkovModel::~HiddenMarkovModel()
	{
		clear();
	}

	HiddenMarkovModel& HiddenMarkovModel::operator = (const HiddenMarkovModel& hmm)
	{
		if (&hmm != this)
		{
			clear();
			copy_(hmm);
			pseudo_counts_ = hmm.pseudo_counts_;
		}
		return *this;
	}

	HMMState* HiddenMarkovModel::getState(const String& name)
	{
		if (name_to_state_.find(name) == name_to_state_.end())
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, name);
		}
		return name_to_state_[name];
	}

	const HMMState* HiddenMarkovModel::getState(const String& name) const
	{
		if (name_to_state_.find(name) == name_to_state_.end())
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, name);
		}
		return name_to_state_.find(name)->second;
	}
	
	void HiddenMarkovModel::writeGraphMLFile(const String& filename)
	{
		set<HMMState*> states = states_;
		map<HMMState*, vector<HMMState*> > transitions;
		
		/*
		stack<HMMState*> s;
		for (map<HMMState*, double>::const_iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
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
				//cerr << ":: " << state->getName() << " " << (*it)->getName() << endl;
				transitions[state].push_back(*it);
				s.push(*it);
			}
		}
		*/

	
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

			set<HMMState*> succ = (*it)->getSuccessorStates();
			for (set<HMMState*>::const_iterator sit = succ.begin(); sit != succ.end(); ++sit)
			{
				transitions[*it].push_back(*sit);
			}
		}


		
		for (map<HMMState*, vector<HMMState*> >::const_iterator it = transitions.begin(); it != transitions.end(); ++it)
		{
			for (vector<HMMState*>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				out << "    <edge source=\"" << it->first->getName() << "\" target=\"" << (*it1)->getName() << "\" directed=\"true\">" << endl;
				out << "      <data key=\"d1\">" << endl;
				out << "        <y:PolyLineEdge>" << endl;
				out << "          <y:EdgeLabel>" << getTransitionProbability(it->first, *it1) << "</y:EdgeLabel>" << endl;
				out << "        </y:PolyLineEdge>" << endl;
				out << "      </data>" << endl;
				out << "    </edge>" << endl;
			}
		}

		out << "  </graph>" << endl;
		out << "</graphml>" << endl;
		
	}

	void HiddenMarkovModel::write(ostream& out) const
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
		for (map<HMMState*, map<HMMState*, double> >::const_iterator it1 = trans_.begin(); it1 != trans_.end(); ++it1)
		{
			for (map<HMMState*, double>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				out << "Transition " << it1->first->getName() << " " << it2->first->getName() << " " << it2->second << endl;
			}
		}
	
		// write synonyms
		/*
		for (map<String, map<String, pair<String, String> > >::const_iterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
		{
			for (map<String, pair<String, String> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				out << "Synonym " << it->first << " " << it2->first << " " << it2->second.first << " " << it2->second.second << endl;
			}
		}*/

    for (map<HMMState*, map<HMMState*, std::pair<HMMState*, HMMState*> > >::const_iterator it1 = synonym_trans_.begin(); it1 != synonym_trans_.end(); ++it1)
    {
      for (map<HMMState*, std::pair<HMMState*, HMMState*> >::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
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
	
	UInt HiddenMarkovModel::getNumberOfStates() const
	{
		return states_.size();
	}

	void HiddenMarkovModel::addNewState(HMMState* s)
	{
		states_.insert(s);
		if (name_to_state_.find(s->getName()) == name_to_state_.end())
		{
			name_to_state_[s->getName()] = s;
		}
		else
		{
			cerr << "HiddenMarkovModel: state name '" << s->getName() << "' (" << s << ") already used!" << endl;
		}
	}

	void HiddenMarkovModel::addNewState(const String& name)
	{
		HMMState* s = new HMMState(name);
		states_.insert(s);
		if (name_to_state_.find(name) == name_to_state_.end())
		{
			name_to_state_[name] = s;
		}
		else
		{
			cerr << "HiddenMarkovModel: state name '" << name << "' (" << name_to_state_[name] << ") already used!" << endl;
		}
		return;
	}

	void HiddenMarkovModel::addSynonymTransition(const String& name1, const String& name2, const String& synonym1, const String& synonym2)
	{
#ifdef SIMPLE_DEBUG2
		cerr << "addSynonymTransition: '" << name1 << "' -> '" << name2 << "' :: '" << synonym1 << "' -> '" << synonym2 << "'" << endl;
#endif
		if (name_to_state_.find(name1) == name_to_state_.end())
		{
			cerr << "state '" << name1 << "' unknown" << endl;
		}
    if (name_to_state_.find(name2) == name_to_state_.end())
    {
      cerr << "state '" << name2 << "' unknown" << endl;
    }
    if (name_to_state_.find(synonym1) == name_to_state_.end())
    {
      cerr << "state '" << synonym1 << "' unknown" << endl;
    }
    if (name_to_state_.find(synonym2) == name_to_state_.end())
    {
      cerr << "state '" << synonym2 << "' unknown" << endl;
    }
		synonym_trans_names_[synonym1][synonym2] = make_pair<String, String>(name1, name2);

		synonym_trans_[name_to_state_[synonym1]][name_to_state_[synonym2]] = make_pair<HMMState*, HMMState*>(name_to_state_[name1], name_to_state_[name2]);

/*
		for (map<String, map<String, std::pair<String, String> > >::const_iterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
    {
      for (map<String, pair<String, String> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
      {
        synonym_trans_[name_to_state_[it->first]][name_to_state_[it2->first]] =
            make_pair<HMMState*, HMMState*>(name_to_state_[it2->second.first], name_to_state_[it2->second.second]);
      }
    }
	*/	
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
#ifdef SIMPLE_DEBUG2
    cerr << "setTransitionProbability: '" << s1 << "' -> '" << s2 << "'" << " " << trans_prob << endl;
#endif

		trans_[name_to_state_[s1]][name_to_state_[s2]] = trans_prob;
		name_to_state_[s1]->addSuccessorState(name_to_state_[s2]);
		name_to_state_[s2]->addPredecessorState(name_to_state_[s1]);
		enabled_trans_[name_to_state_[s1]].insert(name_to_state_[s2]);
		training_steps_count_[name_to_state_[s1]][name_to_state_[s2]] = 0;
	}

	double HiddenMarkovModel::getTransitionProbability(const String& s1, const String& s2) const
	{
		if (name_to_state_.find(s1) == name_to_state_.end())
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, s1);
		}
		if (name_to_state_.find(s2) == name_to_state_.end())
		{
			throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, s2);
		}
		return getTransitionProbability(name_to_state_.find(s1)->second, name_to_state_.find(s2)->second);
	}
	
	double HiddenMarkovModel::getTransitionProbability(HMMState* s1, HMMState* s2) const
	{
		HMMState* state1 = s1;
		HMMState* state2 = s2;

#ifdef SIMPLE_DEBUG2
		cerr << "getTransitionProbability(" << s1->getName() << ", " << s2->getName() << ")" << endl;
#endif
		
		// check if transition is a synonym transition
		if (synonym_trans_.find(state1) != synonym_trans_.end() && 
				synonym_trans_.find(state1)->second.find(state2) != synonym_trans_.find(state1)->second.end())
		{
			HMMState* tmp = synonym_trans_.find(state1)->second.find(state2)->second.first;
			state2 = synonym_trans_.find(state1)->second.find(state2)->second.second; //name_to_state_[p.second];
			state1 = tmp;
		}
		
#ifdef SIMPLE_DEBUG2
		cerr << "getTransitionProbability: " << state1->getName() << " " << state2->getName() << endl;
#endif
		
		// get the transition prob 
		if (trans_.find(state1) != trans_.end() && trans_.find(state1)->second.find(state2) != trans_.find(state1)->second.end())
		{
			return trans_.find(state1)->second.find(state2)->second;
		}
		else
		{
			return 0;
		}
	}

	void HiddenMarkovModel::train()
	{
#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "HiddenMarkovModel::train()" << endl;
#endif
		trained_trans_.clear();

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

		for (map<HMMState*, double>::const_iterator it1 = train_emission_prob_.begin(); it1 != train_emission_prob_.end(); ++it1)
		{
			for (set<HMMState*>::const_iterator it2 = it1->first->getPredecessorStates().begin(); it2 != it1->first->getPredecessorStates().end(); ++it2)
			{
				px += getForwardVariable_(*it2);
			}
		}

		double num_px(0);
		if (px != 0)
		{
			num_px = 1.0/px;
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
			tmp += pseudo_counts_;
			HMMState* s1 = it->first;
			HMMState* s2 = it->second;

			if (synonym_trans_.find(s1) != synonym_trans_.end() && synonym_trans_[s1].find(s2) != synonym_trans_[s1].end())
			{
				HMMState* tmp = synonym_trans_[s1][s2].first;
				s2 = synonym_trans_[s1][s2].second;
				s1 = tmp;
			}
			
			train_count_trans_all_[s1][s2].push_back(tmp);
			if (count_trans_.find(s1) != count_trans_.end() && count_trans_[s1].find(s2) != count_trans_[s1].end())
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
		for (map<HMMState*, double>::iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
		{
			forward_[it->first] = it->second;
		}
		
		for (map<HMMState*, double>::iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
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
		for (map<HMMState*, double>::iterator it = train_emission_prob_.begin(); it != train_emission_prob_.end(); ++it)
		{
			backward_[it->first] = it->second;
		}

		for (map<HMMState*, double>::iterator it = train_emission_prob_.begin(); it != train_emission_prob_.end(); ++it)
		{
			#ifdef HIDDEN_MARKOV_MODEL_DEBUG
			cerr << "b:" << it->first << " " <<  it->first->getName() << " " << it->second << endl;
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
		return forward_.find(state) != forward_.end() ? forward_[state] : 0;
	}

	double HiddenMarkovModel::getBackwardVariable_(HMMState* state)
	{
		return backward_.find(state) != backward_.end() ? backward_[state] : 0;
	}

	void HiddenMarkovModel::evaluate()
	{
		for (map<HMMState*, map<HMMState*, double> >::const_iterator it1 = count_trans_.begin(); it1 != count_trans_.end(); ++it1)
		{
#ifdef EVALUATE_DEBUG
			cerr <<  it1->first->getName() << endl;
#endif
			double sum(0);
			for (map<HMMState*, double>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				if (count_trans_.find(it1->first) != count_trans_.end() && 
						count_trans_[it1->first].find(it2->first) != count_trans_[it1->first].end())
				{
					sum += count_trans_[it1->first][it2->first];
#ifdef EVALUATE_DEBUG
					cerr << it1->first->getName() << " " << it2->first->getName() << " " << count_trans_[it1->first][it2->first] << endl;
#endif
				}
			}

			if (sum != 0)
			{
				for (map<HMMState*, double>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
				{
					if (count_trans_.find(it1->first) != count_trans_.end() && 
							count_trans_[it1->first].find(it2->first) != count_trans_[it1->first].end())
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
#ifdef SIMPLE_DEBUG2
		cerr << "setTrainingEmissionProbability(" << state << "(" << name_to_state_[state] << "), " << prob << ")" << endl;
#endif
		train_emission_prob_[name_to_state_[state]] = prob;
	}

	void HiddenMarkovModel::clearTrainingEmissionProbabilities()
	{
		train_emission_prob_.clear();
	}

	void HiddenMarkovModel::enableTransition(const String& s1, const String& s2)
	{
#ifdef SIMPLE_DEBUG2
		cerr << "enableTransition: '" << s1 << "' -> '" << s2 << "'" << endl;
#endif
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
		for (map<HMMState*, set<HMMState*> >::const_iterator it = enabled_trans_.begin(); it != enabled_trans_.end(); ++it)
		{
			for (set<HMMState*>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				it->first->deleteSuccessorState(*it2);
				(*it2)->deletePredecessorState(it->first);
			}
		}
		enabled_trans_.clear();
	}

	void HiddenMarkovModel::calculateEmissionProbabilities(Map<HMMState*, double>& emission_probs)
	{
		map<HMMState*, double> states = init_prob_;

		while (states.size() != 0)
		{
			map<HMMState*, double> tmp = states;
			for (map<HMMState*, double>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
			{
				for (set<HMMState*>::const_iterator it2 = it->first->getSuccessorStates().begin(); it2 != it->first->getSuccessorStates().end(); ++it2)
				{
					//cerr << "->" << (*it2)->getName() << "=(from " << it->first->getName() << ") " << states[it->first] << " + " << getTransitionProbability(it->first, *it2) << endl;
					if (states.find(*it2) != states.end())
					{
						states[*it2] += states[it->first] * getTransitionProbability(it->first, *it2);
					}
					else
					{
						states[*it2] = states[it->first] * getTransitionProbability(it->first, *it2);
					}
					if (!(*it2)->isHidden())
					{
						if (emission_probs.find(*it2) != emission_probs.end())
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
		for (map<HMMState*, map<HMMState*, double> >::const_iterator it = trans_.begin(); it != trans_.end(); ++it)
		{
			for (map<HMMState*, double>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				cout << it->first->getName() << " -> " << it1->first->getName() << " " << it1->second << " " << training_steps_count_[it->first][it1->first] << ": ";
				vector<double> all_trans = train_count_trans_all_[it->first][it1->first];
				
				if (all_trans.size() != 0)
				{
					double sum = accumulate(all_trans.begin(), all_trans.end(), 0.0);
					double avg(sum/double(all_trans.size()));
					double rsd(0);
					for (UInt i = 0; i != all_trans.size(); ++i)
					{
						cout << all_trans[i] << " ";
						rsd += abs(all_trans[i] - avg);
					}
					cout << "rsd=" << rsd / double(all_trans.size()) / avg;
					cout << ", avg=" << avg;
				}
				
				cout << endl;
			}
		}
		
		cerr << "dump completed" << endl;
	}

	void HiddenMarkovModel::forwardDump()
	{
   set<HMMState*> succ;
    for (map<HMMState*, double>::iterator it = init_prob_.begin(); it != init_prob_.end(); ++it)
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

/*
	void HiddenMarkovModel::buildSynonyms()
	{
		for (map<String, map<String, std::pair<String, String> > >::const_iterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
		{
			for (map<String, pair<String, String> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				synonym_trans_[name_to_state_[it->first]][name_to_state_[it2->first]] = 
						make_pair<HMMState*, HMMState*>(name_to_state_[it2->second.first], name_to_state_[it2->second.second]);
			}
		}
	}*/

	void HiddenMarkovModel::estimateUntrainedTransitions()
	{
		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "Transition training stats" << endl;

		/*
		for (map<HMMState*, map<HMMState*, UInt> >::const_iterator it = training_steps_count_.begin(); it != training_steps_count_.end(); ++it)
		{
			for (map<HMMState*, UInt>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				cerr << it->first->getName() << " " << it1->first->getName() << " " << it1->second << endl;
			}
		}*/

		cerr << "Estimation" << endl;
		#endif


		String residues("ACDEFGHIKMNPQRSTVWY");
		
		vector<String> suffixe;
		HMMState* s2 = name_to_state_["bxyz"];
		HMMState* end_state = name_to_state_["end"];
		for (UInt i = 0; i != residues.size(); ++i)
		{
			s2 = name_to_state_["bxyz"];
			for (UInt j = 0; j != residues.size(); ++j)
			{
				String aa1(residues[i]), aa2(residues[j]);
				if (training_steps_count_[name_to_state_[aa1 + aa2 + "bxyz"]][s2] == 0)
				{
					UInt count(0);
					double sum(0);
					for (UInt k = 0; k != residues.size(); ++k)
					{
						UInt tmp = training_steps_count_[name_to_state_[aa1 + residues[k] + "bxyz"]][s2];
						if (tmp != 0)
						{
							sum += trans_[name_to_state_[aa1 + residues[k] + "bxyz"]][s2];
							count++;
						}
					}
					for (UInt k = 0; k != residues.size(); ++k)
					{
						UInt tmp = training_steps_count_[name_to_state_[residues[k] + aa2 +"bxyz"]][s2];
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

			for (UInt j = 0; j != residues.size(); ++j)
      {
				for (UInt counter = 0; counter < 3; ++counter)
				{
        	String aa1(residues[i]), aa2(residues[j]);
        	if (training_steps_count_[name_to_state_[aa1 + aa2 + "bxyz" + String(counter)]][s2] == 0)
        	{
          	UInt count(0);
          	double sum(0);
          	for (UInt k = 0; k != residues.size(); ++k)
          	{
            	UInt tmp = training_steps_count_[name_to_state_[aa1 + residues[k] + "bxyz" + String(counter)]][s2];
            	if (tmp != 0)
            	{
              	sum += trans_[name_to_state_[aa1 + residues[k] + "bxyz" + String(counter)]][s2];
              	count++;
            	}
          	}
          	for (UInt k = 0; k != residues.size(); ++k)
          	{
            	UInt tmp = training_steps_count_[name_to_state_[residues[k] + aa2 +"bxyz" + String(counter)]][s2];
            	if (tmp != 0)
            	{
              	sum += trans_[name_to_state_[residues[k] + aa2 +"bxyz" + String(counter)]][s2];
              	count++;
            	}
          	}

          	if (count != 0)
          	{
            	#ifdef HIDDEN_MARKOV_MODEL_DEBUG
            	cerr << "setting transitions of " << aa1 << aa2 << "bxyz -> bxyz to " << sum/double(count) << endl;
            	#endif
            	trans_[name_to_state_[aa1 + aa2 + "bxyz" + String(counter)]][s2] = sum/double(count);
            	trans_[name_to_state_[aa1 + aa2 + "bxyz" + String(counter)]][end_state] = 1 - sum/double(count);
          	}
       		}
      	}
			}



			
			
			s2 = name_to_state_["axyz"];
      for (UInt j = 0; j != residues.size(); ++j)
      {
        String aa1(residues[i]), aa2(residues[j]);
        if (training_steps_count_[name_to_state_[aa1 + aa2 + "axyz"]][s2] == 0)
        {
          UInt count(0);
          double sum(0);
          for (UInt k = 0; k != residues.size(); ++k)
          {
            UInt tmp = training_steps_count_[name_to_state_[aa1 + residues[k] + "axyz"]][s2];
            if (tmp != 0)
            {
              sum += trans_[name_to_state_[aa1 + residues[k] + "axyz"]][s2];
              count++;
            }
          }
          for (UInt k = 0; k != residues.size(); ++k)
          {
            UInt tmp = training_steps_count_[name_to_state_[residues[k] + aa2 +"axyz"]][s2];
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
		
			for (UInt j = 0; j != sc_residues.size(); ++j)
			{
				String aa1(residues[i]), sc_res(sc_residues[j]);
				s2 = name_to_state_[sc_res];
				if (training_steps_count_[name_to_state_[aa1 + sc_res]][s2] == 0)
				{
					UInt count(0);
					double sum(0);
					for (UInt k = 0; k != residues.size(); ++k)
					{
						HMMState* s1 = name_to_state_[residues[k] + sc_res];
						UInt tmp = training_steps_count_[s1][s2];
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
        UInt count(0);
        double sum(0);
        for (UInt k = 0; k != residues.size(); ++k)
        {
					HMMState* s1 = name_to_state_[residues[k] + sc_res];
          UInt tmp = training_steps_count_[s1][s2];
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

		s2 = name_to_state_["bk-1"];
		for (UInt i = 0; i != residues.size(); ++i)
		{
			String aa1(residues[i]);
			//cerr << "#training steps " << aa1 << "=" << training_steps_count_[name_to_state_[aa1 + "bk-1"]][s2] << " (" << trans_[name_to_state_[aa1 + "bk-1"]][s2] << ")" << endl;
			if (training_steps_count_[name_to_state_[aa1 + "bk-1"]][s2] == 0)
			{
				UInt count(0);
				double sum(0);
				for (UInt j = 0; j != residues.size(); ++j)
				{
					HMMState* s1 = name_to_state_[residues[j] + "bk-1"];
					UInt tmp = training_steps_count_[s1][s2];
					if (tmp != 0)
					{
						sum += trans_[s1][s2];
						count++;
					}
					//cerr << "Estimating transition of '" << aa1 << "bk-1' -> 'bk-1' to " << sum/(double)count << endl;
					if (count != 0)
					{
						trans_[name_to_state_[aa1 + "bk-1"]][s2] = sum/(double)count;
						trans_[name_to_state_[aa1 + "bk-1"]][end_state] = 1 - sum/(double)count;
					}
				}
			}
		}

		s2 = name_to_state_["bk-2"];
    for (UInt i = 0; i != residues.size(); ++i)
    {
      String aa1(residues[i]);
      if (training_steps_count_[name_to_state_[aa1 + "bk-2"]][s2] == 0)
      {
        UInt count(0);
        double sum(0);
        for (UInt j = 0; j != residues.size(); ++j)
        {
          HMMState* s1 = name_to_state_[residues[j] + "bk-2"];
          UInt tmp = training_steps_count_[s1][s2];
          if (tmp != 0)
          {
            sum += trans_[s1][s2];
            count++;
          }
          if (count != 0)
          {
            trans_[name_to_state_[aa1 + "bk-2"]][s2] = sum/(double)count;
            trans_[name_to_state_[aa1 + "bk-2"]][end_state] = 1 - sum/(double)count;
          }
        }
      }
    }

	}

	void HiddenMarkovModel::setPseudoCounts(double pseudo_counts)
	{
		pseudo_counts_ = pseudo_counts;
		return;
	}

	double HiddenMarkovModel::getPseudoCounts() const
	{
		return pseudo_counts_;
	}

	void HiddenMarkovModel::copy_(const HiddenMarkovModel& source)
	{
		map<HMMState*, HMMState*> old_to_new;
    for (set<HMMState*>::const_iterator it = source.states_.begin(); it != source.states_.end(); ++it)
    {
      HMMState* s = new HMMState(**it);
      states_.insert(s);
      name_to_state_[s->getName()] = s;
      old_to_new[*it] = s;
    }
		
		// trans_
		for (map<HMMState*, map<HMMState*, double> >::const_iterator it1 = source.trans_.begin(); it1 != source.trans_.end(); ++it1)
		{
			for (map<HMMState*, double>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				trans_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
			}
		}

		// count_trans_
    for (map<HMMState*, map<HMMState*, double> >::const_iterator it1 = source.count_trans_.begin(); it1 != source.count_trans_.end(); ++it1)
    {
      for (map<HMMState*, double>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        count_trans_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
      }
    }

		for (map<HMMState*, map<HMMState*, std::vector<double> > >::const_iterator it1 = source.train_count_trans_all_.begin(); it1 != source.train_count_trans_all_.end(); ++it1)
		{
			for (map<HMMState*, vector<double> >::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        train_count_trans_all_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
      }

		}


		for (map<HMMState*, map<HMMState*, UInt> >::const_iterator it1 = source.training_steps_count_.begin(); it1 != source.training_steps_count_.end(); ++it1)
    {
      for (map<HMMState*, UInt>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        training_steps_count_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
      }
    }

		// forward and backward are just temporary objects

		for (map<HMMState*, double>::const_iterator it = source.train_emission_prob_.begin(); it != source.train_emission_prob_.end(); ++it)
		{
			train_emission_prob_[old_to_new[it->first]] = it->second;
		}

		for (map<HMMState*, double>::const_iterator it = source.init_prob_.begin(); it != source.init_prob_.end(); ++it)
		{
			init_prob_[old_to_new[it->first]] = it->second;
		}

		for (std::set<std::pair<HMMState*, HMMState*> >::const_iterator it = source.trained_trans_.begin(); it != source.trained_trans_.end(); ++it)
		{
			trained_trans_.insert(make_pair(old_to_new[it->first], old_to_new[it->second]));
		}

		synonym_trans_names_ = source.synonym_trans_names_;
		pseudo_counts_ = source.pseudo_counts_;


		//buildSynonyms();
    for (map<String, map<String, std::pair<String, String> > >::const_iterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
    {
      for (map<String, pair<String, String> >::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
      {
        synonym_trans_[name_to_state_[it->first]][name_to_state_[it2->first]] =
            make_pair<HMMState*, HMMState*>(name_to_state_[it2->second.first], name_to_state_[it2->second.second]);
      }
    }


		for (map<HMMState*, set<HMMState*> >::const_iterator it1 = source.enabled_trans_.begin(); it1 != source.enabled_trans_.end(); ++it1)
		{
			for (set<HMMState*>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				enabled_trans_[old_to_new[it1->first]].insert(old_to_new[*it2]);
			}
		}

	}
}


