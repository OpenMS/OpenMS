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


#ifndef OPENMS_ANALYSIS_ID_HIDDENMARKOVMODELLIGHT_H
#define OPENMS_ANALYSIS_ID_HIDDENMARKOVMODELLIGHT_H

#include <vector>
#include <set>

#include <OpenMS/DATASTRUCTURES/HashMap.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <utility>

namespace OpenMS 
{
	/** 
	  @brief State for the light weighted Hidden Markov Model implementation

		The lieghweighted Hidden Markov Model implementation has no Strings for
		identifications. Instead of Strings unsigned ints are used for identification
		which can be combined with enums.
	*/
	class HMMStateLight
	{
		public:

			/** @name Constructors and destructors
			 */
			//@{
			/// default constructor with visibility option (default hidden)
			HMMStateLight(bool hidden = true);
			
			/// copy constructor
			HMMStateLight(const HMMStateLight& state);

			/// detailed constructor with id and visibility
			HMMStateLight(Size identifier, bool hidden = true);

			/// destructor
			virtual ~HMMStateLight();
			//@}
		
			///
			void setIdentifier(Size id);
			
			///
			Size getIdentifier() const;

			///
			const HMMStateLight& operator = (const HMMStateLight&);

			///
			void setHidden(bool hidden);

			///
			bool isHidden() const;

			///
			void addPredecessorState(HMMStateLight* state);

			///
			void deletePredecessorState(HMMStateLight* state);

			///
			void addSuccessorState(HMMStateLight* state);

			///
			void deleteSuccessorState(HMMStateLight* state);

			///
			const std::set<HMMStateLight*>& getPredecessorStates() const;

			///
			const std::set<HMMStateLight*>& getSuccessorStates() const;

		protected:
			
			bool hidden_;
			
			Size id_;

			std::set<HMMStateLight*> pre_states_;

			std::set<HMMStateLight*> succ_states_;
	};

				
	class HiddenMarkovModelLight
	{
		/**
		  @brief Hidden Markov Model without Strings which is faster
			
			The light weighted Hidden Markov Model does not make use of Strings. 
			Instead unsigned ints are used which is much faster. 

			Most of the implemented algorithms can only deal with forward connected
			HMMs. The HMM implementation also includes an mechanim for alias states 
			and transitions.
		*/
		public:
			
			/** 
			  @name Constructors and destructors
			 */
			//@{
			/// default constructor
			HiddenMarkovModelLight();

			/// copy constructor
			HiddenMarkovModelLight(const HiddenMarkovModelLight& hmm_new);

			/// destructor
			virtual ~HiddenMarkovModelLight();
			//@}
			
			///
			void writetoYGFFile(const String& filename);
	
			///
			void writeToFile(const String& filename);

			///
			void readFromFile(const String& filename);

			///
			double getTransitionProbability(HMMStateLight*, HMMStateLight*) const;

			///
			double getTransitionProbability(Size id1, Size id2) const;

			///
			void setTransitionProbability(HMMStateLight*, HMMStateLight*, double prob);

			///
			void setTransitionProbability(Size id1, Size id2, double prob);
			
			///
			const Size getNumberOfStates() const;
	
			///
			void addNewState(HMMStateLight*);

			///
			void addSynonymTransition(Size name1, Size name2, Size synonym1, Size synonym2);

			///
			const HiddenMarkovModelLight& operator = (const HiddenMarkovModelLight&);

			///
			void evaluate();

			///
			void train();

			///
			void setTrainingInitialTransitionProbability(HMMStateLight*, double prob);

			///
			void setTrainingInitialTransitionProbability(Size id, double prob);

			///
			void clearTrainingInitialTransitionProbabilities();

			///
			void setTrainingEmissionProbability(HMMStateLight*, double prob);

			///
			void setTrainingEmissionProbability(Size, double prob);

			///
			void clearTrainingEmissionProbabilities();

			///
			void enableTransition(HMMStateLight*, HMMStateLight*);

			///
			void enableTransition(Size id1, Size id2);

			///
			void disableTransition(HMMStateLight*, HMMStateLight*);

			///
			void disableTransition(Size id1, Size id2);

			///
			void disableTransitions();

			///
			void calculateEmissionProbabilities(HashMap<HMMStateLight*, double>& emission_probs);

			///
			void dump();

			///
			void forwardDump(); 

			///
			void buildSynonyms();

			///
			void estimateUntrainedTransitions();

			///
			HMMStateLight* getState(Size id1);

			///
			const HMMStateLight* getState(Size id1) const;

			///
			void addIdToName(Size id, const String& name);

		protected:
			
			void calculateForwardPart_();

			void calculateBackwardPart_();

			double getForwardVariable_(HMMStateLight*);

			double getBackwardVariable_(HMMStateLight*);

		private:

			HashMap<HMMStateLight*, HashMap<HMMStateLight*, double> > trans_;

			HashMap<HMMStateLight*, HashMap<HMMStateLight*, double> > count_trans_;

			HashMap<HMMStateLight*, HashMap<HMMStateLight*, double> > train_count_trans_;

			HashMap<HMMStateLight*, HashMap<HMMStateLight*, Size> > training_steps_count_;

			HashMap<HMMStateLight*, double> forward_;

			HashMap<HMMStateLight*, double> backward_;

			HashMap<Size, HMMStateLight*> id_to_state_;

			HashMap<HMMStateLight*, double> train_emission_prob_;

			HashMap<HMMStateLight*, double> init_train_prob_;

			std::set<HMMStateLight*> states_;

			std::set<std::pair<HMMStateLight*, HMMStateLight*> > trained_trans_;

			HashMap<Size, HashMap<Size, std::pair<Size, Size> > > synonym_trans_names_;

			HashMap<HMMStateLight*, HashMap<HMMStateLight*, std::pair<HMMStateLight*, HMMStateLight*> > > synonym_trans_;

			HashMap<HMMStateLight*, std::set<HMMStateLight*> > enabled_trans_;

			HashMap<Size, String> id_to_name_;
	};
}
#endif

