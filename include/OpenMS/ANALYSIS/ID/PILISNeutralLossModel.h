// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_ID_PILISNEUTRALLOSSMODEL_H
#define OPENMS_ANALYSIS_ID_PILISNEUTRALLOSSMODEL_H

#include <vector>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS 
{
	class AASequence;
	
	/** 
	  @brief This class implements the simulation of the spectra from PILIS

		PILIS uses a HMM based structure to model the population of fragment ions
		from a peptide. The spectrum generator can be accessed via the getSpectrum
		method.
		 
		@htmlinclude OpenMS_PILISNeutralLossModel.parameters

		@ingroup Analysis_ID
	*/	
	class OPENMS_DLLAPI PILISNeutralLossModel : public DefaultParamHandler
	{
		friend class PILISNeutralLossModelGenerator;

		public:
			
			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor
			PILISNeutralLossModel();

			/// copy constructor
			PILISNeutralLossModel(const PILISNeutralLossModel& model);
			
			/// destructor
			virtual ~PILISNeutralLossModel();
			//@}

			/// assignment operator
			PILISNeutralLossModel& operator = (const PILISNeutralLossModel& mode);
			
			/** @name Accessors
			*/
			//@{
			/// performs a training step; needs as parameters a spectrum with annotated sequence and charge; returns the intensity sum of the matched peaks
			DoubleReal train(const RichPeakSpectrum& spec, const AASequence& peptide, DoubleReal ion_weight, UInt charge, DoubleReal peptide_weight);

			/// given a peptide (a ion) the model returns the peaks with intensities relative to initial_prob
			void getIons(std::vector<RichPeak1D>& peaks, const AASequence& peptide, DoubleReal initial_prob);

			/// sets the hidden markov model
			void setHMM(const HiddenMarkovModel& model);
			
			/// writes the HMM to the given file in the GraphML format. A detailed description of the GraphML format can be found under http://graphml.graphdrawing.org/
			const HiddenMarkovModel& getHMM() const;
			
			/// generates the models
			void generateModel();
			
			/// this method evaluates the model after training; it should be called after all training steps with train
			void evaluate();
			//@}

		protected:

			/// extracts the precursor and related intensities of a training spectrum
			DoubleReal getIntensitiesFromSpectrum_(const RichPeakSpectrum& train_spec, Map<String, DoubleReal>& pre_ints, DoubleReal ion_weight, const AASequence& peptide, UInt charge);

			/// trains precursor and related peaks
			void trainIons_(DoubleReal initial_probability, const Map<String, DoubleReal>& intensities, const AASequence& peptide);
			
			/// estimates the precursor intensities 
			void getIons_(Map<String, DoubleReal>& intensities, DoubleReal initial_probability, const AASequence& precursor);

			/// enables the states needed for precursor training/simulation
			void enableIonStates_(const AASequence& peptide);

			/// precursor model used
			HiddenMarkovModel hmm_precursor_;

			/// 
			UInt num_explicit_;

			void updateMembers_();
	};
}
#endif

