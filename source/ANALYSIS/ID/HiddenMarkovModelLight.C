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


#include <OpenMS/ANALYSIS/ID/HiddenMarkovModelLight.h>

#include <iostream>
#include <fstream>
#include <stack>
#include <cstdlib>

//#define PSEUDO_COUNTS 0.0000001
#define PSEUDO_COUNTS   1e-15

#define HIDDEN_MARKOV_MODEL_DEBUG
#undef HIDDEN_MARKOV_MODEL_DEBUG

using namespace std;

namespace OpenMS 
{
	// HMMSTATE 
	HMMStateLight::HMMStateLight()
		: hidden_(true)
	{
	}

	HMMStateLight::HMMStateLight(UInt id, bool hidden)
		: hidden_(hidden),
			id_(id)
	{
	}

	HMMStateLight::HMMStateLight(const HMMStateLight& state)
		: hidden_(state.hidden_),
			id_(state.id_)
	{
	}
	
	HMMStateLight::~HMMStateLight()
	{
	}
	
	void HMMStateLight::setIdentifier(UInt id)
	{
		id_ = id;
	}

	UInt HMMStateLight::getIdentifier() const
	{
		return id_;
	}

	HMMStateLight& HMMStateLight::operator = (const HMMStateLight& state)
	{
		hidden_ = state.hidden_;
		id_ = state.id_;
		pre_states_.clear();
		succ_states_.clear();
		//pre_states_ = state.pre_states_;
		//succ_states_ = state.succ_states_;
		return *this;
	}

	void HMMStateLight::addSuccessorState(HMMStateLight* state)
	{
		succ_states_.insert(state);
	}

	void HMMStateLight::deleteSuccessorState(HMMStateLight* state)
	{
		succ_states_.erase(state);
	}

	void HMMStateLight::addPredecessorState(HMMStateLight* state)
	{
		pre_states_.insert(state);
	}

	void HMMStateLight::deletePredecessorState(HMMStateLight* state)
	{
		pre_states_.erase(state);
	}

	const set<HMMStateLight::HMMStateLight*>& HMMStateLight::getPredecessorStates() const
	{
		return pre_states_;
	}

	const set<HMMStateLight::HMMStateLight*>& HMMStateLight::getSuccessorStates() const
	{
		return succ_states_;
	}

	bool HMMStateLight::isHidden() const
	{
		return hidden_;
	}

	void HMMStateLight::setHidden(bool hidden)
	{
		hidden_ = hidden;
	}

	// The hidden markov model
	HiddenMarkovModelLight::HiddenMarkovModelLight()
		:	pseudo_counts_(PSEUDO_COUNTS)
	{
	}

	HiddenMarkovModelLight::HiddenMarkovModelLight(const HiddenMarkovModelLight& hmm)
/*
 * :	trans_(hmm.trans_),
      count_trans_(hmm.count_trans_),
      train_count_trans_(hmm.train_count_trans_),
      training_steps_count_(hmm.training_steps_count_),
      forward_(hmm.forward_),
      backward_(hmm.backward_),
      id_to_state_(hmm.id_to_state_),
      train_emission_prob_(hmm.train_emission_prob_),
      init_train_prob_(hmm.init_train_prob_),
      trained_trans_(hmm.trained_trans_),
      synonym_trans_names_(hmm.synonym_trans_names_),
      synonym_trans_(hmm.synonym_trans_),
      enabled_trans_(hmm.enabled_trans_),
			pseudo_counts_(hmm.pseudo_counts_)*/
	{
	/*
		for (set<HMMStateLight*>::const_iterator it = hmm.states_.begin(); it != hmm.states_.end(); ++it)
		{
			HMMStateLight* s = new HMMStateLight(**it);
			states_.insert(s);
			id_to_state_[s->getIdentifier()] = s;
		}
		*/
		copy_(hmm);
	}
	
	HiddenMarkovModelLight::~HiddenMarkovModelLight()
	{
		clear();
		/*
		for (set<HMMStateLight*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
		{
			delete *it;
		}
		*/
	}

	HiddenMarkovModelLight& HiddenMarkovModelLight::operator = (const HiddenMarkovModelLight& hmm)
	{
		if (&hmm != this)
		{
			clear();
			copy_(hmm);
	
		/*
			for (set<HMMStateLight*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
			{
				delete *it;
			}

			for (set<HMMStateLight*>::const_iterator it = hmm.states_.begin(); it != hmm.states_.end(); ++it)
			{
				HMMStateLight* s = new HMMStateLight(**it);
				states_.insert(s);
				id_to_state_[s->getIdentifier()] = s;
			}

			trans_ = hmm.trans_;
			count_trans_ = hmm.count_trans_;
			train_count_trans_ = hmm.train_count_trans_;
			training_steps_count_ = hmm.training_steps_count_;
			forward_ = hmm.forward_;
			backward_ = hmm.backward_;
			id_to_state_ = hmm.id_to_state_;
			train_emission_prob_ = hmm.train_emission_prob_;
			init_train_prob_ = hmm.init_train_prob_;
			trained_trans_ = hmm.trained_trans_;
			synonym_trans_names_ = hmm.synonym_trans_names_;
			synonym_trans_ = hmm.synonym_trans_;
			enabled_trans_ = hmm.enabled_trans_;
			pseudo_counts_ = hmm.pseudo_counts_;
			*/
		}
		return *this;
	}

	HMMStateLight* HiddenMarkovModelLight::getState(UInt id)
	{
		return id_to_state_[id];
	}

	const HMMStateLight* HiddenMarkovModelLight::getState(UInt id) const
	{
		return id_to_state_[id];
	}
	
	void HiddenMarkovModelLight::writeGraphMLFile(const String& filename)
	{
		set<HMMStateLight*> states;
		Map<HMMStateLight*, vector<HMMStateLight*> > transitions;
		
		stack<HMMStateLight*> s;
		for (Map<HMMStateLight*, double>::ConstIterator it = init_train_prob_.begin(); it != init_train_prob_.end(); ++it)
		{
			s.push(it->first);
		}

		while (s.size() != 0)
		{
			HMMStateLight* state = s.top();
			states.insert(state);
			s.pop();
			for (set<HMMStateLight*>::const_iterator it = state->getSuccessorStates().begin(); it != state->getSuccessorStates().end(); ++it)
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
		for (set<HMMStateLight*>::const_iterator it = states.begin(); it != states.end(); ++it)
		{
			out << "    <node id=\"" << id_to_name_[(*it)->getIdentifier()] << "\">" << endl;
			out << "      <data key=\"d0\">" << endl;
			out << "        <y:ShapeNode>" << endl;
			out << "          <y:NodeLabel>" << id_to_name_[(*it)->getIdentifier()] << "</y:NodeLabel>" << endl;
			out << "        </y:ShapeNode>" << endl;
			out << "      </data>" << endl;
			out << "    </node>" << endl;
		}

		for (Map<HMMStateLight*, vector<HMMStateLight*> >::ConstIterator it = transitions.begin(); it != transitions.end(); ++it)
		{
			for (vector<HMMStateLight*>::const_iterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				out << "    <edge source=\"" << id_to_name_[it->first->getIdentifier()] << "\" target=\"" << id_to_name_[(*it1)->getIdentifier()] << "\">" << endl;
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

	void HiddenMarkovModelLight::readFromFile(const String& filename)
	{
		// clear old model
		for (set<HMMStateLight*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
		{
			delete *it;
		}
		trans_.clear();
		count_trans_.clear();
		forward_.clear();
		backward_.clear();
		id_to_state_.clear();
		train_emission_prob_.clear();
		init_train_prob_.clear();
		states_.clear();
		trained_trans_.clear();
		synonym_trans_.clear();
		synonym_trans_names_.clear();
		
		ifstream in(filename.c_str());
		char buf[10000];
		HMMStateLight* s = 0;
		Map<HMMStateLight*, vector<UInt> > pre_states;
		Map<HMMStateLight*, vector<UInt> > succ_states;
		while (in.getline(buf, 10000, '\n'))
		{
			//cerr << ">> " << buf << endl;
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
				s = new HMMStateLight((UInt)split[1].toInt());
			}
			if (split_0 == "hidden_:")
			{
				if (split[1] == "false")
				{
					s->setHidden(false);
				}
			
				states_.insert(s);
				id_to_state_[s->getIdentifier()] = s;
			}
		}

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "states block ready (#states=" << states_.size() << ", id_to_state_.size()=" << id_to_state_.size() << ")" << endl;
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
			trans_[id_to_state_[split[0].toInt()]][id_to_state_[split[2].toInt()]] = atof(split[3].c_str());
		}

		// synonym trans block
		while (in.getline(buf, 10000, '\n'))
		{
			String line(buf);
			line.trim();
			vector<String> split;
			line.split(' ', split);
			synonym_trans_names_[split[0].toInt()][split[1].toInt()] = make_pair<UInt, UInt>(split[2].toInt(), split[3].toInt());
		}

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		cerr << "trans_.size()=" << trans_.size() << ", synonym_trans_names_.size()=" << synonym_trans_names_.size() << endl;
		#endif

		buildSynonyms();

		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		dump();
		#endif
	}

  void HiddenMarkovModelLight::clear()
  {
    for (set<HMMStateLight*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
    {
      delete *it;
    }
    trans_.clear();
    count_trans_.clear();
    forward_.clear();
    backward_.clear();
    id_to_state_.clear();
    train_emission_prob_.clear();
    init_train_prob_.clear();
    states_.clear();
    trained_trans_.clear();
    synonym_trans_.clear();
    synonym_trans_names_.clear();
    return;
  }


	UInt HiddenMarkovModelLight::getNumberOfStates() const
	{
		return states_.size();
	}

	void HiddenMarkovModelLight::addNewState(HMMStateLight* s)
	{
		states_.insert(s);
		if (!id_to_state_.has(s->getIdentifier()))
		{
			id_to_state_[s->getIdentifier()] = s;
		}
		else
		{
			cerr << "HiddenMarkovModelLight: state name '" << s->getIdentifier() << "' (" << s << ") already used!" << endl;
		}
	}

	void HiddenMarkovModelLight::addSynonymTransition(UInt name1, UInt name2, UInt synonym1, UInt synonym2)
	{
		synonym_trans_names_[synonym1][synonym2] = make_pair<UInt, UInt>(name1, name2);
	}

	void HiddenMarkovModelLight::setTransitionProbability(HMMStateLight * s1, HMMStateLight * s2, double trans_prob)
	{
		trans_[s1][s2] = trans_prob;
		s1->addSuccessorState(s2);
		s2->addPredecessorState(s1);
		enabled_trans_[s1].insert(s2);
		training_steps_count_[s1][s2] = 0;
	}

	void HiddenMarkovModelLight::setTransitionProbability(UInt s1, UInt s2, double trans_prob)
	{
		trans_[id_to_state_[s1]][id_to_state_[s2]] = trans_prob;
		id_to_state_[s1]->addSuccessorState(id_to_state_[s2]);
		id_to_state_[s2]->addPredecessorState(id_to_state_[s1]);
		enabled_trans_[id_to_state_[s1]].insert(id_to_state_[s2]);
		training_steps_count_[id_to_state_[s1]][id_to_state_[s2]] = 0;
	}

  double HiddenMarkovModelLight::getTransitionProbability(UInt s1, UInt s2) const
  {
    if (!id_to_state_.has(s1))
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,String(s1));
    }
    if (!id_to_state_.has(s2))
    {
      throw Exception::ElementNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__,String(s2));
    }
    return getTransitionProbability(id_to_state_[s1], id_to_state_[s2]);
  }

	double HiddenMarkovModelLight::getTransitionProbability(HMMStateLight* s1, HMMStateLight* s2) const
	{
		HMMStateLight* state1 = s1;
		HMMStateLight* state2 = s2;
		
		// check if transition is a synonym transition
		if (synonym_trans_.has(state1) && synonym_trans_[state1].has(state2))
		{
			HMMStateLight* tmp = synonym_trans_[state1][state2].first;
			state2 = synonym_trans_[state1][state2].second; //id_to_state_[p.second];
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

	void HiddenMarkovModelLight::train()
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

		for (Map<HMMStateLight*, double>::ConstIterator it1 = train_emission_prob_.begin(); it1 != train_emission_prob_.end(); ++it1)
		{
			for (set<HMMStateLight*>::const_iterator it2 = it1->first->getPredecessorStates().begin(); it2 != it1->first->getPredecessorStates().end(); ++it2)
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
		for (set<pair<HMMStateLight*, HMMStateLight*> >::const_iterator it = trained_trans_.begin(); it != trained_trans_.end(); ++it)
		{
			double tmp(0);
			tmp = num_px * getForwardVariable_(it->first) * getBackwardVariable_(it->second) * getTransitionProbability(it->first, it->second);
			tmp += pseudo_counts_/*PSEUDO_COUNTS*/;
			HMMStateLight* s1 = it->first;
			HMMStateLight* s2 = it->second;
			train_count_trans_[s1][s2] = tmp;

			if (synonym_trans_.has(s1) && synonym_trans_[s1].has(s2))
			{
				HMMStateLight* tmp = synonym_trans_[s1][s2].first; //id_to_state_[p.first];
				s2 = synonym_trans_[s1][s2].second; // id_to_state_[p.second];
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

	void HiddenMarkovModelLight::calculateForwardPart_()
	{
		forward_.clear();
		set<HMMStateLight*> succ;
		for (Map<HMMStateLight*, double>::Iterator it = init_train_prob_.begin(); it != init_train_prob_.end(); ++it)
		{
			forward_[it->first] = it->second;
		}
		
		for (Map<HMMStateLight*, double>::Iterator it = init_train_prob_.begin(); it != init_train_prob_.end(); ++it)
		{
			succ.insert(it->first->getSuccessorStates().begin(), it->first->getSuccessorStates().end());
		
			while (succ.size() != 0)
			{
				set<HMMStateLight*> succ_new;
				for (set<HMMStateLight*>::const_iterator it = succ.begin(); it != succ.end(); ++it)
				{
					set<HMMStateLight*> pre = (*it)->getPredecessorStates();
					double sum(0);
					for (set<HMMStateLight*>::const_iterator it2 = pre.begin(); it2 != pre.end(); ++it2)
					{
						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
						cerr << "\tadding from " << (*it2)->getIdentifier() << " " << getForwardVariable_(*it2) << " * " << getTransitionProbability(*it2, *it) << endl;
						#endif
						sum += getForwardVariable_(*it2) * getTransitionProbability(*it2, *it);
						trained_trans_.insert(make_pair<HMMStateLight*, HMMStateLight*>(*it2, *it));
					}
					forward_[*it] = sum;
					succ_new.insert((*it)->getSuccessorStates().begin(), (*it)->getSuccessorStates().end());
					#ifdef HIDDEN_MARKOV_MODEL_DEBUG
					cerr << "f: " << (*it)->getIdentifier() << "\t" << sum << endl;
					#endif
				}
				succ = succ_new;
			}
		}
		#ifdef HIDDEN_MARKOV_MODEL_DEBUG
		for (Map<HMMStateLight*, double>::ConstIterator it = forward_.begin(); it != forward_.end(); ++it)
		{
			cerr << it->first->getIdentifier() << "\t" << it->second << endl;
		}
		#endif
	}

	void HiddenMarkovModelLight::calculateBackwardPart_()
	{
		backward_.clear();
		set<HMMStateLight*> pre;
		for (Map<HMMStateLight*, double>::Iterator it = train_emission_prob_.begin(); it != train_emission_prob_.end(); ++it)
		{
			backward_[it->first] = it->second;
		}

		for (Map<HMMStateLight*, double>::Iterator it = train_emission_prob_.begin(); it != train_emission_prob_.end(); ++it)
		{
			#ifdef HIDDEN_MARKOV_MODEL_DEBUG
			cerr << "b:" << it->first->getIdentifier() << " " << it->second << endl;
			#endif
			pre.insert(it->first->getPredecessorStates().begin(), it->first->getPredecessorStates().end());
		
			while (pre.size() != 0)
			{
				set <HMMStateLight*> pre_new;
				for (set<HMMStateLight*>::const_iterator it = pre.begin(); it != pre.end(); ++it)
				{
					set<HMMStateLight*> succ = (*it)->getSuccessorStates();
					double sum(0);
					for (set<HMMStateLight*>::const_iterator it2 = succ.begin(); it2 != succ.end(); ++it2)
					{
						#ifdef HIDDEN_MARKOV_MODEL_DEBUG
						cerr << "\tadding from " << (*it2)->getIdentifier() << " " << getBackwardVariable_(*it2) << " * " << getTransitionProbability(*it, *it2) << endl;
						#endif
						sum += getBackwardVariable_(*it2) * getTransitionProbability(*it, *it2);
						trained_trans_.insert(make_pair<HMMStateLight*, HMMStateLight*>(*it, *it2));
					}
					backward_[*it] = sum;
					pre_new.insert((*it)->getPredecessorStates().begin(), (*it)->getPredecessorStates().end());
					#ifdef HIDDEN_MARKOV_MODEL_DEBUG
					cerr << "b: " << (*it)->getIdentifier() << "\t" << sum << endl;
					#endif
				}
				pre = pre_new;
			}
		}
	}

	double HiddenMarkovModelLight::getForwardVariable_(HMMStateLight* state)
	{
		return forward_.has(state) ? forward_[state] : 0;
	}

	double HiddenMarkovModelLight::getBackwardVariable_(HMMStateLight* state)
	{
		return backward_.has(state) ? backward_[state] : 0;
	}

	void HiddenMarkovModelLight::evaluate()
	{
		for (Map<HMMStateLight*, Map<HMMStateLight*, double> >::ConstIterator it1 = count_trans_.begin(); it1 != count_trans_.end(); ++it1)
		{
			double sum(0);
			for (Map<HMMStateLight*, double>::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				if (count_trans_.has(it1->first) && count_trans_[it1->first].has(it2->first))
				{
					sum += count_trans_[it1->first][it2->first];
				}
			}

			if (sum != 0)
			{
				for (Map<HMMStateLight*, double>::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
				{
					if (count_trans_.has(it1->first) && count_trans_[it1->first].has(it2->first))
					{
						trans_[it1->first][it2->first] = count_trans_[it1->first][it2->first] / sum;
					}
				}
			}
		}
	}

	void HiddenMarkovModelLight::setInitialTransitionProbability(UInt state, double prob)
	{
		init_train_prob_[id_to_state_[state]] = prob;
	}

	void HiddenMarkovModelLight::clearInitialTransitionProbabilities()
	{
		init_train_prob_.clear();
	}

	void HiddenMarkovModelLight::setTrainingEmissionProbability(UInt state, double prob)
	{
		train_emission_prob_[id_to_state_[state]] = prob;
	}

	void HiddenMarkovModelLight::clearTrainingEmissionProbabilities()
	{
		train_emission_prob_.clear();
	}

	void HiddenMarkovModelLight::enableTransition(UInt s1, UInt s2)
	{
		enableTransition(id_to_state_[s1], id_to_state_[s2]);
	}

	void HiddenMarkovModelLight::enableTransition(HMMStateLight* s1, HMMStateLight* s2)
	{
		s1->addSuccessorState(s2);
		s2->addPredecessorState(s1);
		enabled_trans_[s1].insert(s2);
	}

	void HiddenMarkovModelLight::disableTransition(UInt s1, UInt s2)
	{
		disableTransition(id_to_state_[s1], id_to_state_[s2]);
	}
	
	void HiddenMarkovModelLight::disableTransition(HMMStateLight* s1, HMMStateLight* s2)
	{
		s1->deleteSuccessorState(s2);
		s2->deletePredecessorState(s1);
		enabled_trans_[s1].erase(s2);
	}

	void HiddenMarkovModelLight::disableTransitions()
	{
		for (Map<HMMStateLight*, set<HMMStateLight*> >::ConstIterator it = enabled_trans_.begin(); it != enabled_trans_.end(); ++it)
		{
			for (set<HMMStateLight*>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				//disableTransition(it->first, *it2);
				it->first->deleteSuccessorState(*it2);
				(*it2)->deletePredecessorState(it->first);
			}
		}
		enabled_trans_.clear();
	}

	void HiddenMarkovModelLight::calculateEmissionProbabilities(Map<HMMStateLight*, double>& emission_probs)
	{
		Map<HMMStateLight*, double> states = init_train_prob_;

		while (states.size() != 0)
		{
			Map<HMMStateLight*, double> tmp = states;
			for (Map<HMMStateLight*, double>::ConstIterator it = tmp.begin(); it != tmp.end(); ++it)
			{
				for (set<HMMStateLight*>::const_iterator it2 = it->first->getSuccessorStates().begin(); it2 != it->first->getSuccessorStates().end(); ++it2)
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

	void HiddenMarkovModelLight::dump()
	{
		cout << "dump of transitions: " << endl;
		for (Map<HMMStateLight*, Map<HMMStateLight*, double> >::ConstIterator it = trans_.begin(); it != trans_.end(); ++it)
		{
			for (Map<HMMStateLight*, double>::ConstIterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				cout << it->first->getIdentifier() << " -> " << it1->first->getIdentifier() << " " << it1->second << endl;
			}
		}
		cout << "dump completed" << endl;
	}

	void HiddenMarkovModelLight::forwardDump()
	{
   set<HMMStateLight*> succ;
    for (Map<HMMStateLight*, double>::Iterator it = init_train_prob_.begin(); it != init_train_prob_.end(); ++it)
    {
      succ.insert(it->first->getSuccessorStates().begin(), it->first->getSuccessorStates().end());

      while (succ.size() != 0)
      {
        set<HMMStateLight*> succ_new;
        for (set<HMMStateLight*>::const_iterator it = succ.begin(); it != succ.end(); ++it)
        {
					cerr << (*it)->getIdentifier() << endl;
          succ_new.insert((*it)->getSuccessorStates().begin(), (*it)->getSuccessorStates().end());
        }
        succ = succ_new;
      }
    }
	}

	void HiddenMarkovModelLight::buildSynonyms()
	{
		for (Map<UInt, Map<UInt, std::pair<UInt,  UInt> > >::ConstIterator it = synonym_trans_names_.begin(); it != synonym_trans_names_.end(); ++it)
		{
			for (Map<UInt, pair<UInt, UInt> >::ConstIterator it2 = it->second.begin(); it2 != it->second.end(); ++it2)
			{
				synonym_trans_[id_to_state_[it->first]][id_to_state_[it2->first]] = 
						make_pair<HMMStateLight*, HMMStateLight*>(id_to_state_[it2->second.first], id_to_state_[it2->second.second]);
			}
		}
	}

	void HiddenMarkovModelLight::estimateUntrainedTransitions()
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		/*
		for (Map<HMMStateLight*, Map<HMMStateLight*, UInt> >::ConstIterator it = training_steps_count_.begin(); it != training_steps_count_.end(); ++it)
		{
			for (Map<HMMStateLight*, UInt>::ConstIterator it1 = it->second.begin(); it1 != it->second.end(); ++it1)
			{
				cerr << it->first->getIdentifier() << " " << it1->first->getIdentifier() << " " << it1->second << endl;
			}
		}

		cerr << "Estimation" << endl;
		String residues("ACDEFGHIKMNPQRSTVWY");
		
		vector<String> suffixe;
		HMMStateLight* s2 = id_to_state_["bxyz"];
		HMMStateLight* end_state = id_to_state_["end"];
		for (UInt i = 0; i != residues.size(); ++i)
		{
			s2 = id_to_state_["bxyz"];
			for (UInt j = 0; j != residues.size(); ++j)
			{
				String aa1(residues[i]), aa2(residues[j]);
				if (training_steps_count_[id_to_state_[aa1 + aa2 + "bxyz"]][s2] == 0)
				{
					UInt count(0);
					double sum(0);
					for (UInt k = 0; k != residues.size(); ++k)
					{
						UInt tmp = training_steps_count_[id_to_state_[aa1 + residues[k] + "bxyz"]][s2];
						if (tmp != 0)
						{
							sum += trans_[id_to_state_[aa1 + residues[k] + "bxyz"]][s2];
							count++;
						}
					}
					for (UInt k = 0; k != residues.size(); ++k)
					{
						UInt tmp = training_steps_count_[id_to_state_[residues[k] + aa2 +"bxyz"]][s2];
						if (tmp != 0)
						{
							sum += trans_[id_to_state_[residues[k] + aa2 +"bxyz"]][s2];
							count++;
						}
					}

					if (count != 0)
					{
						cerr << "setting transitions of " << aa1 << aa2 << "bxyz -> bxyz to " << sum/double(count) << endl;
						trans_[id_to_state_[aa1 + aa2 + "bxyz"]][s2] = sum/double(count);
						trans_[id_to_state_[aa1 + aa2 + "bxyz"]][end_state] = 1 - sum/double(count);
					}
				}
			}

			s2 = id_to_state_["axyz"];
      for (UInt j = 0; j != residues.size(); ++j)
      {
        String aa1(residues[i]), aa2(residues[j]);
        if (training_steps_count_[id_to_state_[aa1 + aa2 + "axyz"]][s2] == 0)
        {
          UInt count(0);
          double sum(0);
          for (UInt k = 0; k != residues.size(); ++k)
          {
            UInt tmp = training_steps_count_[id_to_state_[aa1 + residues[k] + "axyz"]][s2];
            if (tmp != 0)
            {
              sum += trans_[id_to_state_[aa1 + residues[k] + "axyz"]][s2];
              count++;
            }
          }
          for (UInt k = 0; k != residues.size(); ++k)
          {
            UInt tmp = training_steps_count_[id_to_state_[residues[k] + aa2 +"axyz"]][s2];
            if (tmp != 0)
            {
              sum += trans_[id_to_state_[residues[k] + aa2 +"axyz"]][s2];
              count++;
            }
          }
					if (count != 0)
					{
	          cerr << "setting transitions of " << aa1 << aa2 << "axyz -> axyz to " << sum/double(count) << endl;
  	        trans_[id_to_state_[aa1 + aa2 + "axyz"]][s2] = sum/double(count);
    	      trans_[id_to_state_[aa1 + aa2 + "axyz"]][end_state] = 1 - sum/double(count);
					}
        }
      }

			// sc and cr
			String sc_residues("HKDE");
		
			for (UInt j = 0; j != sc_residues.size(); ++j)
			{
				String aa1(residues[i]), sc_res(sc_residues[j]);
				s2 = id_to_state_[sc_res];
				if (training_steps_count_[id_to_state_[aa1 + sc_res]][s2] == 0)
				{
					UInt count(0);
					double sum(0);
					for (UInt k = 0; k != residues.size(); ++k)
					{
						HMMStateLight* s1 = id_to_state_[residues[k] + sc_res];
						UInt tmp = training_steps_count_[s1][s2];
						if (tmp != 0)
						{
							sum += trans_[s1][s2];
							count++;
						}
					}

					if (count != 0)
					{
						trans_[id_to_state_[aa1 + sc_res]][s2] = sum/double(count);
						trans_[id_to_state_[aa1 + sc_res]][end_state] = 1 - sum/double(count);
					}
				}
			}

			String aa1(residues[i]), sc_res("RSC");
			s2 = id_to_state_["R"];

      if (training_steps_count_[id_to_state_[aa1 + sc_res]][s2] == 0)
      {
        UInt count(0);
        double sum(0);
        for (UInt k = 0; k != residues.size(); ++k)
        {
					HMMStateLight* s1 = id_to_state_[residues[k] + sc_res];
          UInt tmp = training_steps_count_[s1][s2];
          if (tmp != 0)
          {
            sum += trans_[s1][s2];
            count++;
          }
        }
				if (count != 0)
				{
	        trans_[id_to_state_[aa1 + sc_res]][s2] = sum/double(count);
	      	trans_[id_to_state_[aa1 + sc_res]][end_state] = 1 - sum/double(count);
				}
			}
		}
		*/
	}

	void HiddenMarkovModelLight::addIdToName(UInt id, const String& name)
	{
		id_to_name_[id] = name;
	}

	void HiddenMarkovModelLight::write(ostream& out)
	{
		//ofstream outfile(filename.c_str());
		// states
		for (set<HMMStateLight*>::const_iterator it = states_.begin(); it != states_.end(); ++it)
		{
			out << "State " << id_to_name_[(*it)->getIdentifier()];
			if (!(*it)->isHidden())
			{
				out << " false";
			}
			out << endl;
		}
		// transitions
		for (Map<HMMStateLight*, Map<HMMStateLight*, double> >::ConstIterator it1 = trans_.begin(); it1 != trans_.end(); ++it1)
		{
			for (Map<HMMStateLight*, double>::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				out << "Transition " << id_to_name_[it1->first->getIdentifier()] << " " << id_to_name_[it2->first->getIdentifier()] << " " << it2->second << endl;
			}
		}
		// synonym transitions
		for (Map<HMMStateLight*, Map<HMMStateLight*, std::pair<HMMStateLight*, HMMStateLight*> > >::ConstIterator it1 = synonym_trans_.begin(); it1 != synonym_trans_.end(); ++it1)
		{
			for (Map<HMMStateLight*, std::pair<HMMStateLight*, HMMStateLight*> >::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
			{
				out << "Synonym " << id_to_name_[it1->first->getIdentifier()] << " " << id_to_name_[it2->first->getIdentifier()] << " " << id_to_name_[it2->second.first->getIdentifier()] << " " << id_to_name_[it2->second.second->getIdentifier()] << endl;
			}
		}
	}
	
	void HiddenMarkovModelLight::setPseudoCounts(double pseudo_counts)
	{
		pseudo_counts_ = pseudo_counts;
	}

	double HiddenMarkovModelLight::getPseudoCounts() const
	{
		return pseudo_counts_;
	}

  void HiddenMarkovModelLight::copy_(const HiddenMarkovModelLight& source)
  {
    Map<HMMStateLight*, HMMStateLight*> old_to_new;
    for (set<HMMStateLight*>::const_iterator it = source.states_.begin(); it != source.states_.end(); ++it)
    {
      HMMStateLight* s = new HMMStateLight(**it);
      states_.insert(s);
      id_to_state_[s->getIdentifier()] = s;
      old_to_new[*it] = s;
    }

    // trans_
    for (Map<HMMStateLight*, Map<HMMStateLight*, double> >::ConstIterator it1 = source.trans_.begin(); it1 != source.trans_.end(); ++it1)
    {
      for (Map<HMMStateLight*, double>::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        trans_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
      }
    }

    // count_trans_
    for (Map<HMMStateLight*, Map<HMMStateLight*, double> >::ConstIterator it1 = source.count_trans_.begin(); it1 != source.count_trans_.end(); ++it1)
    {
      for (Map<HMMStateLight*, double>::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        count_trans_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
      }
    }

    for (Map<HMMStateLight*, Map<HMMStateLight*, UInt> >::ConstIterator it1 = source.training_steps_count_.begin(); it1 != source.training_steps_count_.end(); ++it1)
    {
      for (Map<HMMStateLight*, UInt>::ConstIterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        training_steps_count_[old_to_new[it1->first]][old_to_new[it2->first]] = it2->second;
      }
    }

		// forward and backward are just temporary objects

    for (Map<HMMStateLight*, double>::ConstIterator it = source.train_emission_prob_.begin(); it != source.train_emission_prob_.end(); ++it)
    {
      train_emission_prob_[old_to_new[it->first]] = it->second;
    }

    for (Map<HMMStateLight*, double>::ConstIterator it = source.init_train_prob_.begin(); it != source.init_train_prob_.end(); ++it)
    {
      init_train_prob_[old_to_new[it->first]] = it->second;
    }

    for (std::set<std::pair<HMMStateLight*, HMMStateLight*> >::const_iterator it = source.trained_trans_.begin(); it != source.trained_trans_.end(); ++it)
    {
      trained_trans_.insert(make_pair(old_to_new[it->first], old_to_new[it->second]));
    }

    synonym_trans_names_ = source.synonym_trans_names_;
    pseudo_counts_ = source.pseudo_counts_;

    buildSynonyms();

    for (Map<HMMStateLight*, set<HMMStateLight*> >::ConstIterator it1 = source.enabled_trans_.begin(); it1 != source.enabled_trans_.end(); ++it1)
    {
      for (set<HMMStateLight*>::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2)
      {
        enabled_trans_[old_to_new[it1->first]].insert(old_to_new[*it2]);
      }
    }

		id_to_name_ = source.id_to_name_;
	}

}

