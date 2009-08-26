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
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_ID_HIDDENMARKOVMODEL_H
#define OPENMS_ANALYSIS_ID_HIDDENMARKOVMODEL_H

#include <vector>
#include <set>

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <utility>

namespace OpenMS 
{
	/** 
	  @brief Hidden Markov Model State class for the Hidden Markov Model
	*/
	class OPENMS_DLLAPI HMMState
	{
		public:

			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor
			HMMState();
			
			/// copy constructor
			HMMState(const HMMState& state);

			/// constructor with name and visibility option
			HMMState(const String& name, bool hidden = true);

			/// destructor
			virtual ~HMMState();
			//@}
		
			///
			HMMState& operator = (const HMMState&);
		
			/** Accessors
			*/
			//@{
			/// sets the name of the state
			void setName(const String& name);
			
			/// returns the name of the state
			const String& getName() const;

			/// sets the hidden property to the state
			void setHidden(bool hidden);

			/// returns true if the state is hidden
			bool isHidden() const;

			/// adds the given predecessor state to the list
			void addPredecessorState(HMMState* state);

			/// deletes the given predecessor state from the list
			void deletePredecessorState(HMMState* state);

			/// add the given successor state to the list
			void addSuccessorState(HMMState* state);

			/// deletes the given successor state from the list
			void deleteSuccessorState(HMMState* state);

			/// returns the predecessor states of the state
			const std::set<HMMState*>& getPredecessorStates() const;

			/// return the successor states of the state
			const std::set<HMMState*>& getSuccessorStates() const;
			//@}

		protected:
	
			///
			bool hidden_;
			
			///
			String name_;

			///
			std::set<HMMState*> pre_states_;

			///
			std::set<HMMState*> succ_states_;
	};

	
	/** 
	  @brief 	Hidden Markov Model implementation of PILIS
		
						Hidden Markov Model implementation suitable for forward conncected HMMs.
						The HMM is mostly used within PILIS. For further details have a look at
						the docs of PILIS.
	*/
	class OPENMS_DLLAPI HiddenMarkovModel
	{
		public:
					
			/** @name Constructors and destructors
			 */
			//@{
			/// default constructor
			HiddenMarkovModel();

			/// copy constructor
			HiddenMarkovModel(const HiddenMarkovModel& hmm_new);

			/// destructor
			virtual ~HiddenMarkovModel();
			//@}

			///
			HiddenMarkovModel& operator = (const HiddenMarkovModel&);
			
			/** Accessors
			*/
			//@{
			/** @brief writes the HMM into a file in GraphML format

					A detailed description of the GraphML format can be found under
					http://graphml.graphdrawing.org/
			*/
			void writeGraphMLFile(const String& filename);
	
			/// writes the HMM into an outstream
			void write(std::ostream& out) const;

			/// read a HMM from the given file
			//void readFromFile(const String& filename);

			/// returns the transition probability of the given states
			double getTransitionProbability(HMMState* s1, HMMState* s2) const;

			/// returns the transition probability of the given state names
			double getTransitionProbability(const String& s1, const String& s2) const;

			/// sets the transition probability of the given states to prob
			void setTransitionProbability(HMMState* s1, HMMState* s2, double prob);

			/// sets the transition probability of the given state names to prob
			void setTransitionProbability(const String& s1, const String& s2, double prob);
			
			/// return the number of states
			Size getNumberOfStates() const;

			/// registers a new state to the HMM
			void addNewState(HMMState* state);

			/// registers a new state to the HMM
			void addNewState(const String& name);

			/// add a new synonym transition to the given state names
			void addSynonymTransition(const String& name1, const String& name2, const String& synonym1, const String& synonym2);

			/// evaluate the HMM, estimates the transition probabilities from the training
			void evaluate();

			/// trains the HMM; initial probabilities and emission probabilities of the emitting states should be set
			void train();

			/// sets the initial transition probability of the given state to prob
			void setInitialTransitionProbability(const String& state, double prob);

			/// clears the initial probabilities 
			void clearInitialTransitionProbabilities();

			/// sets the emission probability of the given state to prob
			void setTrainingEmissionProbability(HMMState* state, double prob);

			/// sets the emission probability of the given state to prob
			void setTrainingEmissionProbability(const String& state, double prob);

			/// clear the emission probabilities
			void clearTrainingEmissionProbabilities();

			/// enables a transition; adds s1 to the predecessor list of s2 and s2 to the successor list of s1
			void enableTransition(HMMState* s1, HMMState* s2);

			/// enables a transition; adds s1 to the predecessor list of s2 and s2 to the successor list of s1
			void enableTransition(const String& s1, const String& s2);

			/// disables the transition; deletes the nodes from the predeccessor/successor list repsectively
			void disableTransition(HMMState* s1, HMMState* s2);

			/// disables the transition; deletes the nodes from the predeccessor/successor list repsectively
			void disableTransition(const String& s1, const String& s2);

			/// disables all transitions
			void disableTransitions();

			/// calculates the emission probabilities of the HMM (of course only of the non-hidden states)
			void calculateEmissionProbabilities(Map<HMMState*, double>& emission_probs);

			/// writes some stats to cerr
			void dump();

			/// writes some info of the forward "matrix" to cerr
			void forwardDump(); 

			/// builds a synonyms structure, needed when synonyms are used
			//void buildSynonyms();

			/// estimates the transition probabilities of not trained transitions by averages similar trained ones
			void estimateUntrainedTransitions();

			/// returns the state with the given name
			HMMState* getState(const String& name);

			/// returns the state with the given name
			const HMMState* getState(const String& name) const;

			/// clears all data
			void clear();

			/// 
			void setPseudoCounts(double pseudo_counts);

			///
			double getPseudoCounts() const;

			void setVariableModifications(const StringList& modifications);
			//@}
			
		protected:
			
			/// performs the forward algorithm
			void calculateForwardPart_();

			/// performs the backward algorithm
			void calculateBackwardPart_();

			/// returns the forward variable 
			double getForwardVariable_(HMMState*);

			/// returns the backward variable
			double getBackwardVariable_(HMMState*);

		private:

			// transition probs
			Map<HMMState*, Map<HMMState*, double> > trans_;

			// transition prob counts
			Map<HMMState*, Map<HMMState*, double> > count_trans_;

			Map<HMMState*, Map<HMMState*, std::vector<double> > > count_trans_all_;

			// all transition probs of all training steps (for model checking)
			Map<HMMState*, Map<HMMState*, std::vector<double> > > train_count_trans_all_;

			// number of training steps of the transitions
			Map<HMMState*, Map<HMMState*, Size> > training_steps_count_;

			// forward variables
			Map<HMMState*, double> forward_;

			// backward variables
			Map<HMMState*, double> backward_;

			// name to state Mapping
			Map<String, HMMState*> name_to_state_;

			// emission probabilities
			Map<HMMState*, double> train_emission_prob_;

			// initial transition probabilities
			Map<HMMState*, double> init_prob_;

			// all states of the HMM
			std::set<HMMState*> states_;

			// trained transitions
			std::set<std::pair<HMMState*, HMMState*> > trained_trans_;

			// synonym transitions Mapping
			Map<String, Map<String, std::pair<String, String> > > synonym_trans_names_;

			// synonym transitions
			Map<HMMState*, Map<HMMState*, std::pair<HMMState*, HMMState*> > > synonym_trans_;

			// transitions which are enabled
			Map<HMMState*, std::set<HMMState*> > enabled_trans_;

			// pseudocounts used in this instance
			double pseudo_counts_;

			// copy all the stuff from one HMM to this 
			void copy_(const HiddenMarkovModel& source);

			StringList var_modifications_;
	};
}
#endif

