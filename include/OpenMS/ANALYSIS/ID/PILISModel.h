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
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

namespace OpenMS 
{
	class PILISModel 
	{
		public:
						
			PILISModel();

			PILISModel(const PILISModel& model);
			
			virtual ~PILISModel();

			void initModel();

			// training
			void train(const PeakSpectrum&, const AASequence& peptide, UnsignedInt charge);

			void readFromFile(const String& file_name);

			void writetoYGFFile(const String& filename);

			void writeToFile(const String& file_name);

			// TODO
			void getSpectrumAlignment(HashMap<Size, Size>& peak_map, const PeakSpectrum& spec1, const PeakSpectrum& spec2);

			void getSpectrum(PeakSpectrum& spec, const AASequence& peptide, UnsignedInt charge);

			void evaluate();


		protected:
		
			enum NeutralLossType_
			{
				LOSS_TYPE_H2O = 0,
				LOSS_TYPE_NH3,
				LOSS_TYPE_NONE
			};

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
			
			struct PrecursorPeaks_
			{
				double pre;
				double pre_H2O;
				double pre_NH3;
			};
		
			struct IonPeaks_
			{
				HashMap<IonType_, std::vector<double> > ints;
			};
			
			
			void initModel_();
	
			void initPrecursorModel_();

			void initLossModels_();

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

			
			// 
			void getPrecursorIntensitiesFromSpectrum_(const PeakSpectrum& train_spec, PrecursorPeaks_& pre_ints, double peptide_weight, UnsignedInt charge);

			double getIntensitiesFromSpectrum_(const PeakSpectrum& train_spec, IonPeaks_& ion_ints, const AASequence& peptide, UnsignedInt charge);
			
			// internal train methods
			void trainPrecursorIons_(double initial_probability, double intensity, double intensity_NH3, double intensity_H2O, const AASequence& peptide);

			void trainNeutralLossesFromIon_(double initial_probability, const HashMap<NeutralLossType_, double>& intensities, Residue::ResidueType ion_type, double ion_intensity, const AASequence& ion);


			// internal getter methods
			//double getPrecursorPeakIntensity_(double initial_probability, const AASequence& peptide);

			void getPrecursorIons_(HashMap<NeutralLossType_, double>& intensities, double initial_probability, const AASequence& precursor);

			void getNeutralLossesFromIon_(HashMap<NeutralLossType_, double>& intensities, double initial_probability, Residue::ResidueType ion_type, const AASequence& ion);


			// enables the states needed
			void enablePrecursorIonStates_(const AASequence& peptide);

			void enableNeutralLossStates_(Residue::ResidueType ion_type, const AASequence& ion);

			
			// get the initial transition probabilities from the proton dist, returns true if charge remote is enabled
			bool getInitialTransitionProbabilities_(std::vector<double>& bb_init, 
																							std::vector<double>& cr_init, 
																							std::vector<double>& sc_init, 
																							const HashMap<Size, double>& bb_charges,
																							const HashMap<Size, double>& sc_charges,
																							const AASequence& peptide);

			// add peaks to spectrum
			void addPeaks_(double mz, int charge, double mz_offset, double intensity, PeakSpectrum& spectrum, const IsotopeDistribution& id, const String& name);
			
			// parse model file
			void parseModelFile(const String& filename, HiddenMarkovModelLight* model);

			static ResidueDB res_db_;

			HiddenMarkovModel hmm_;

			HiddenMarkovModelLight hmm_precursor_;

			HashMap<Residue::ResidueType, HiddenMarkovModelLight> hmms_losses_;

			ProtonDistributionModel prot_dist_;

			TheoreticalSpectrumGenerator tsg_;

			HashMap<String, States_> name_to_enum_;

			HashMap<States_, String> enum_to_name_;

	};
}
#endif

