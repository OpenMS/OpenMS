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


#ifndef OPENMS_ANALYSIS_ID_HIDDENMARKOVMODEL_H
#define OPENMS_ANALYSIS_ID_HIDDENMARKOVMODEL_H

#include <vector>
#include <set>

#include <OpenMS/DATASTRUCTURES/HashMap.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <utility>

namespace OpenMS 
{
	/** 
	  @brief Hidden Markov Model State class for the Hidden Markov Model
	*/
	class HMMState
	{
		public:

			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor with option to create hidden/visible state
			HMMState(bool hidden = true);
			
			/// copy constructor
			HMMState(const HMMState& state);

			/// constructor with name and visibility option
			HMMState(const String& name, bool hidden = true);

			/// destructor
			virtual ~HMMState();
			//@}
			
			///
			void setName(const String&);
			
			///
			const String& getName() const;

			/// 
			const HMMState& operator = (const HMMState&);

			///
			double getEmissionProbability() const;

			///
			void setEmissionProbability(double ep);

			///
			void setHidden(bool hidden);

			///
			bool isHidden() const;

			///
			void addPredecessorState(HMMState* state);

			///
			void deletePredecessorState(HMMState* state);

			///
			void addSuccessorState(HMMState* state);

			///
			void deleteSuccessorState(HMMState* state);

			///
			const std::set<HMMState*>& getPredecessorStates() const;

			///
			const std::set<HMMState*>& getSuccessorStates() const;

		protected:
	
			///
			bool hidden_;
			
			///
			String name_;

			///
			double emission_prob_;

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
	class HiddenMarkovModel
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
			void writetoYGFFile(const String& filename);
	
			///
			void writeToFile(const String& filename);

			///
			void readFromFile(const String& filename);

			///
			double getTransitionProbability(HMMState*, HMMState*) const;

			///
			double getTransitionProbability(const String&, const String&) const;

			///
			void setTransitionProbability(HMMState*, HMMState*, double prob);

			///
			void setTransitionProbability(const String&, const String&, double prob);
			
			///
			const Size getNumberOfStates() const;

			///
			void addNewState(HMMState*);

			//void addSynonymState(const String& name, const String& synonym);

			///
			void addSynonymTransition(const String& name1, const String& name2, const String& synonym1, const String& synonym2);

			///
			const HiddenMarkovModel& operator = (const HiddenMarkovModel&);

			///
			void evaluate();

			///
			void train();

			///
			void setTrainingInitialTransitionProbability(HMMState*, double prob);

			///
			void setTrainingInitialTransitionProbability(const String&, double prob);

			///
			void clearTrainingInitialTransitionProbabilities();

			///
			void setTrainingEmissionProbability(HMMState*, double prob);

			///
			void setTrainingEmissionProbability(const String&, double prob);

			///
			void clearTrainingEmissionProbabilities();

			///
			void enableTransition(HMMState*, HMMState*);

			///
			void enableTransition(const String&, const String&);

			///
			void disableTransition(HMMState*, HMMState*);

			///
			void disableTransition(const String&, const String&);

			///
			void disableTransitions();

			///
			void calculateEmissionProbabilities(HashMap<HMMState*, double>& emission_probs);

			///
			void dump();

			///
			void forwardDump(); 

			///
			void buildSynonyms();

			///
			void estimateUntrainedTransitions();

			///
			HMMState* getState(const String& name);

			///
			const HMMState* getState(const String& name) const;
			
		protected:
			
			///
			void calculateForwardPart_();

			///
			void calculateBackwardPart_();

			///
			double getForwardVariable_(HMMState*);

			///
			double getBackwardVariable_(HMMState*);

		private:

			HashMap<HMMState*, HashMap<HMMState*, double> > trans_;

			HashMap<HMMState*, HashMap<HMMState*, double> > count_trans_;

			HashMap<HMMState*, HashMap<HMMState*, double> > train_count_trans_;

			HashMap<HMMState*, HashMap<HMMState*, Size> > training_steps_count_;

			HashMap<HMMState*, double> forward_;

			HashMap<HMMState*, double> backward_;

			HashMap<String, HMMState*> name_to_state_;

			HashMap<HMMState*, double> train_emission_prob_;

			HashMap<HMMState*, double> init_train_prob_;

			std::set<HMMState*> states_;

			std::set<std::pair<HMMState*, HMMState*> > trained_trans_;

			HashMap<String, HashMap<String, std::pair<String, String> > > synonym_trans_names_;

			HashMap<HMMState*, HashMap<HMMState*, std::pair<HMMState*, HMMState*> > > synonym_trans_;

			HashMap<HMMState*, std::set<HMMState*> > enabled_trans_;
	};
}
#endif

