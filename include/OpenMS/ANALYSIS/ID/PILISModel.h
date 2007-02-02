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


#ifndef OPENMS_ANALYSIS_ID_PILISMODEL_H
#define OPENMS_ANALYSIS_ID_PILISMODEL_H

#include <OpenMS/CHEMISTRY/Residue.h>
#include <vector>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/HashMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/ANALYSIS/ID/HiddenMarkovModel.h>
#include <OpenMS/ANALYSIS/ID/HiddenMarkovModelLight.h>
#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS 
{
	class AASequence;
	
	/** 
	  @brief This class implements the simulation of the spectra from PILIS

		PILIS uses a HMM based structure to model the population of fragment ions
		from a peptide. The spectrum generator can be accessed via the getSpectrum
		method.
	*/	
	class PILISModel : public DefaultParamHandler
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
			void train(const PeakSpectrum&, const AASequence& peptide, UnsignedInt charge);

			/** reads the model parameters from the given files
			    @param filename filename of the model
			*/ 
			void readFromFile(const String& filename);

			/// writes the HMM to the given file in the YGF format
			void writetoYGFFile(const String& filename);

			/** writes the model parameters into the given files
			    @param filename filename of the base model
			*/			 
			void writeToFile(const String& filename);

			/// greedy specturm aligner, should be replaced by a better algorithm
			//void getSpectrumAlignment(HashMap<Size, Size>& peak_map, const PeakSpectrum& spec1, const PeakSpectrum& spec2);

			/// simulates a spectrum with the model of the given peptide and charge and writes it to the given PeakSpectrum
			void getSpectrum(PeakSpectrum& spec, const AASequence& peptide, UnsignedInt charge);

			/// this method evaluates the model after training; it should be called after all training steps with train
			void evaluate();

			/// mutable access to the parameters
			//Param& getParam();

			/// non-mutable access to the parameters
			//const Param& getParam() const;

			/// set the parameters
			//void setParam(const Param& param);

			/// resets the parameters to the default values
			//void resetToDefaultParam();
			//@}

		protected:

			/// enumerates the basic loss types implemented
			enum NeutralLossType_
			{
				LOSS_TYPE_H2O = 0,
				LOSS_TYPE_NH3,
				LOSS_TYPE_NONE
			};

			/// enumeration of the basic ion types used
			enum IonType_
			{
				AIon = 0,
				BIon,
				YIon,
				BIon_H2O,
				BIon_NH3,
				YIon_H2O,
				YIon_NH3
			};
			
			/// describes precursor and related peaks
			struct PrecursorPeaks_
			{
				double pre;
				double pre_H2O;
				double pre_NH3;
			};
		
			/// describes ions peaks and the relatives of them
			struct IonPeaks_
			{
				HashMap<IonType_, std::vector<double> > ints;
			};
			
			/// initializes the model
			void initModels_();
	
			/// the states of the precursor and loss states
			enum States_
			{
				PRE_MH = 0,
				PRE_MH_H2O,
				PRE_MH_NH3,
				PRE_END,
				PRE_ION,
				PRE_BASE1,
				PRE_BASE2,
				PRE_H2O_S,
				PRE_H2O_T,
				PRE_H2O_E,
				PRE_H2O_D,
				PRE_H2O_Q1,
				PRE_H2O_CTERM,
				PRE_NH3_K,
				PRE_NH3_R,
				B_H2O,
				B_NH3,
				B_LOSS_END,
				B_ION,
				B_BASE1,
				B_BASE2,
				B_H2O_S,
				B_H2O_T,
				B_H2O_E,
				B_H2O_D,
				B_NH3_K,
				B_NH3_R,
        Y_H2O,
        Y_NH3,
        Y_LOSS_END,
        Y_ION,
        Y_BASE1,
        Y_BASE2,
        Y_H2O_S,
        Y_H2O_T,
        Y_H2O_E,
        Y_H2O_D,
				Y_H2O_Q1,
				Y_H2O_CTERM,
        Y_NH3_K,
        Y_NH3_R
      };			

			
			/// extracts the precursor and related intensities of a training spectrum
			void getPrecursorIntensitiesFromSpectrum_(const PeakSpectrum& train_spec, PrecursorPeaks_& pre_ints, double peptide_weight, UnsignedInt charge);

			/// extracts the ions intensities of a training spectrum
			double getIntensitiesFromSpectrum_(const PeakSpectrum& train_spec, IonPeaks_& ion_ints, const AASequence& peptide, UnsignedInt charge);

			/// aligns two spectra a writes the intensities from the first which matches the second to the vector
			double getIntensitiesFromComparison_(const PeakSpectrum& train_spec, const PeakSpectrum& theo_spec, std::vector<double>& intensities);

			/// trains precursor and related peaks
			void trainPrecursorIons_(double initial_probability, double intensity, double intensity_NH3, double intensity_H2O, const AASequence& peptide);

			/// trains neutral losses an related peaks
			void trainNeutralLossesFromIon_(double initial_probability, const HashMap<NeutralLossType_, double>& intensities, Residue::ResidueType ion_type, double ion_intensity, const AASequence& ion);

			/// estimates the precursor intensities 
			void getPrecursorIons_(HashMap<NeutralLossType_, double>& intensities, double initial_probability, const AASequence& precursor);

			/// estimates the neutral losses of an ion
			void getNeutralLossesFromIon_(HashMap<NeutralLossType_, double>& intensities, double initial_probability, Residue::ResidueType ion_type, const AASequence& ion);

			/// enables the states needed for precursor training/simulation
			void enablePrecursorIonStates_(const AASequence& peptide);

			/// enables the states needed for neutral loss training/simulation
			void enableNeutralLossStates_(Residue::ResidueType ion_type, const AASequence& ion);

			/// get the initial transition probabilities from the proton dist, returns true if charge remote is enabled
			bool getInitialTransitionProbabilities_(std::vector<double>& bb_init, 
																							std::vector<double>& cr_init, 
																							std::vector<double>& sc_init, 
																							const HashMap<Size, double>& bb_charges,
																							const HashMap<Size, double>& sc_charges,
																							const AASequence& peptide);

			/// add peaks to spectrum
			void addPeaks_(double mz, int charge, double mz_offset, double intensity, PeakSpectrum& spectrum, const IsotopeDistribution& id, const String& name);
		
			/// parse the base model
			void parseBaseModel_(const TextFile::ConstIterator& begin, const TextFile::ConstIterator& end);
			
			/// parse model file of losses and precursor models
			void parseHMMLightModel_(const TextFile::ConstIterator& begin, const TextFile::ConstIterator& end, HiddenMarkovModelLight& model);

			/// residue db used
			static ResidueDB res_db_;

			/// base model used
			HiddenMarkovModel hmm_;

			/// precursor model used
			HiddenMarkovModelLight hmm_precursor_;

			/// loss models used
			HashMap<Residue::ResidueType, HiddenMarkovModelLight> hmms_losses_;

			/// proton distribution model
			ProtonDistributionModel prot_dist_;

			/// theoretical spectrum generator (needed for training/aligning and spectrum intensity extraction)
			TheoreticalSpectrumGenerator tsg_;

			/// name to enum mapping of the losses/precursor states
			HashMap<String, States_> name_to_enum_;

			/// enum to name mapping of the losses/precursor states
			HashMap<States_, String> enum_to_name_;

			/// params_ of the model
			//Param param_;

			/// default parameters of the model
			//Param default_;

			/// true if the instance is valid
			bool valid_;

			/// stores the peaks of a spectrum
			HashMap<double, std::vector<Peak> > peaks_;

			/// the alignment algorithm used
			SpectrumAlignment spectra_aligner_;
	};
}
#endif

