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
			double train(const RichPeakSpectrum&, const AASequence& peptide, double ion_weight, UInt charge, double peptide_weight);

			/** reads the model parameters from the given files
			    @param filename filename of the model
			*/ 
			//void readFromFile(const String& filename);

			/// 
			void getIons(std::vector<RichPeak1D>& peaks, const AASequence& peptide, double initial_prob);

			///
			void setHMM(const HiddenMarkovModel& model);
			/// writes the HMM to the given file in the GraphML format. A detailed description of the GraphML format can be found under http://graphml.graphdrawing.org/
			const HiddenMarkovModel& getHMM() const;
			//void writeGraphMLFile(const String& filename);

			/** writes the model parameters into the given files
			    @param filename filename of the base model
			*/			 
			//void writeToFile(const String& filename);
			void generateModel();
			
			
			/// this method evaluates the model after training; it should be called after all training steps with train
			void evaluate();
			//@}

		protected:

			/// extracts the precursor and related intensities of a training spectrum
			double getIntensitiesFromSpectrum_(const RichPeakSpectrum& train_spec, Map<String, double>& pre_ints, double ion_weight, const AASequence& peptide, UInt charge);

			/// trains precursor and related peaks
			void trainIons_(double initial_probability, const Map<String, double>& intensities, const AASequence& peptide);
			
			/// estimates the precursor intensities 
			void getIons_(Map<String, double>& intensities, double initial_probability, const AASequence& precursor);

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

