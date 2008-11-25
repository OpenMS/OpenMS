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


#ifndef OPENMS_ANALYSIS_ID_PILISMODEL_H
#define OPENMS_ANALYSIS_ID_PILISMODEL_H

#include <vector>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>
#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/ID/PILISNeutralLossModel.h>


namespace OpenMS 
{
	class AASequence;
	
	/** 
	  @brief This class implements the simulation of the spectra from PILIS

		PILIS uses a HMM based structure to model the population of fragment ions
		from a peptide. The spectrum generator can be accessed via the getSpectrum
		method.
		 
		@htmlinclude OpenMS_PILISModel.parameters

		@ingroup Analysis_ID
	*/	
	class PILISModel
		: public DefaultParamHandler
	{
		public:
						
			/** @name Constructors and destructors
			*/
			//@{
			/// default constructor
			PILISModel();

			/// copy constructor
			PILISModel(const PILISModel& model);
			
			/// destructor
			virtual ~PILISModel();
			//@}

			/// assignment operator
			PILISModel& operator = (const PILISModel& mode);
			
			/** @name Accessors
			*/
			//@{
			/// performs a training step; needs as parameters a spectrum with annotated sequence and charge
			void train(const RichPeakSpectrum&, const AASequence& peptide, UInt charge);

			/** reads the model parameters from the given files
			    @param filename filename of the model
			*/ 
			void readFromFile(const String& filename);

			/// writes the HMM to the given file in the GraphML format. A detailed description of the GraphML format can be found under http://graphml.graphdrawing.org/
			void writeGraphMLFile(const String& filename);

			/** writes the model parameters into the given files
			    @param filename filename of the base model
			*/			 
			void writeToFile(const String& filename);


			// 
			void init(bool generate_models = true);
			
			/// greedy specturm aligner, should be replaced by a better algorithm
			//void getSpectrumAlignment(Map<UInt, UInt>& peak_map, const PeakSpectrum& spec1, const PeakSpectrum& spec2);

			/// simulates a spectrum with the model of the given peptide and charge and writes it to the given PeakSpectrum
			void getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, UInt charge);

			/// sets the main hidden Markov model
			void setHMM(const HiddenMarkovModel& model);

			/// set the precursor hidden Markov model
			void setPrecursorHMM(const HiddenMarkovModel& model);
			
			/// this method evaluates the model after training; it should be called after all training steps with train
			void evaluate();
			//@}

		protected:

			/// enumeration of the basic ion types used
			enum IonType_
			{
				AIon = 0,
				BIon,
				B2Ion,
				YIon,
				BIon_H2O,
				BIon_NH3,
				YIon_H2O,
				YIon_NH3
			};
			
			/// describes ions peaks and the relatives of them
			struct IonPeaks_
			{
				Map<IonType_, std::vector<double> > ints;
			};
			
			/// extracts the ions intensities of a training spectrum
			double getIntensitiesFromSpectrum_(const RichPeakSpectrum& train_spec, IonPeaks_& ion_ints, const AASequence& peptide, UInt charge);

			/// aligns two spectra a writes the intensities from the first which matches the second to the vector
			double getIntensitiesFromComparison_(const RichPeakSpectrum& train_spec, const RichPeakSpectrum& theo_spec, std::vector<double>& intensities);

			/// trains neutral losses an related peaks
			void trainNeutralLossesFromIon_(double initial_probability, const Map<String, double>& intensities, IonType_ ion_type, double ion_intensity, const AASequence& ion);

			/// estimates the neutral losses of an ion
			void getNeutralLossesFromIon_(Map<String, double>& intensities, double initial_probability, IonType_ ion_type, const AASequence& ion);

			/// enables the states needed for neutral loss training/simulation
			void enableNeutralLossStates_(IonType_ ion_type, const AASequence& ion);

			/// get the initial transition probabilities from the proton dist, returns true if charge remote is enabled
			bool getInitialTransitionProbabilities_(std::vector<double>& bb_init, 
																							std::vector<double>& cr_init, 
																							std::vector<double>& sc_init, 
																							const Map<UInt, double>& bb_charges,
																							const Map<UInt, double>& sc_charges,
																							const AASequence& peptide);

			double getAvailableBackboneCharge_(const AASequence& ion, Residue::ResidueType res_type, int charge);

			/// add peaks to spectrum
			void addPeaks_(double mz, int charge, double mz_offset, double intensity, RichPeakSpectrum& spectrum, const IsotopeDistribution& id, const String& name);
		
			/// parse the base model
			void parseHMMModel_(const TextFile::ConstIterator& begin, const TextFile::ConstIterator& end, HiddenMarkovModel& hmm);
		
			/// base model used
			HiddenMarkovModel hmm_;

			/// proton distribution model
			ProtonDistributionModel prot_dist_;

			/// theoretical spectrum generator (needed for training/aligning and spectrum intensity extraction)
			TheoreticalSpectrumGenerator tsg_;

			/// true if the instance is valid
			bool valid_;

			/// stores the peaks of a spectrum
			Map<double, std::vector<RichPeak1D> > peaks_;

			/// the alignment algorithm used
			SpectrumAlignment spectra_aligner_;

			/// precursor model used
			PILISNeutralLossModel precursor_model_cr_;

			PILISNeutralLossModel precursor_model_cd_;

			PILISNeutralLossModel a_ion_losses_cr_;
			PILISNeutralLossModel a_ion_losses_cd_;
			
			PILISNeutralLossModel b_ion_losses_cr_;
			PILISNeutralLossModel b_ion_losses_cd_;

			PILISNeutralLossModel b2_ion_losses_cr_;
			PILISNeutralLossModel b2_ion_losses_cd_;
			
			PILISNeutralLossModel y_ion_losses_cr_;
			PILISNeutralLossModel y_ion_losses_cd_;

			void updateMembers_();
	};
}
#endif

