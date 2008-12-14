// -*- mode: C++; tab-width: 2; -*-
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


#ifndef OPENMS_ANALYSIS_ID_PROTONDISTRIBUTIONMODEL_H
#define OPENMS_ANALYSIS_ID_PROTONDISTRIBUTIONMODEL_H

#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <vector>

namespace OpenMS 
{
	class AASequence;

	/**
		@brief A proton distribution model to calculate the proton distribution over charged peptides
 	
		The model uses proton affinity values of backbone nitrogens and sidechains to calculate the 
		proton distribution of charged peptide among these sites. The possible sites are the peptide
		bonds between the amino acids, the side chains and the C-terminus and N-terminus. The calculation
		is done calculating a Boltzmann distribution of the sites.

		Details and the proton affinities can be found in
		Z. Zhang, Prediction of Low-Energy Collision-Induced Dissociation Spectra of Peptides,
    Anal. Chem., 76 (14), 3908 - 3922, 2004

		A proton distribution can be calculated using the getProtonDistribution method. The backbone 
		probabilities are reported in the first parameter (index 0 for the N-terminus, index 1 for the 
		first peptide bond...), the site chain probabilities are reported in the second parameter 
		(index 0, for the first amino acid...). The peptide and the number of protons as well as type 
		of peptide (can be Reside::YIon for peptides and y-ions and any other ion type).

		Charge state intensities of differently charged equal (e.g. y7+ and y7++) ions can be calculated
		using the getChargeStateIntensities function.
	 
		@htmlinclude OpenMS_ProtonDistributionModel.parameters
	*/
	class OPENMS_DLLAPI ProtonDistributionModel : public DefaultParamHandler
	{
		public:
			
			/** Constructor and destructors
			*/
			//@{
			/// default constructor
			ProtonDistributionModel();

			/// copy constructor
			ProtonDistributionModel(const ProtonDistributionModel& model);
			
			/// destructor
			virtual ~ProtonDistributionModel();
			//@}

			/// assignment operator 
			ProtonDistributionModel& operator = (const ProtonDistributionModel& pdm);

			/** @name Enumerations
			*/
			//@{
			/// the type of fragmentation
      enum FragmentationType
      {
        ChargeDirected = 0,
        ChargeRemote,
        SideChain
      };
			//@}

			/** @brief calculates a proton distribution of the given charged peptide
			
					@param bb_charges the calculated probabilities of the backbone sites (including N-terminus and C-terminus)
					@param sc_charges the calculated probabilities of the side chain sites
					@param peptide the peptide
					@param charge the charge
					@param res_type the type of the ion given in peptide. Peptides are handled as y-ions, i.e. Residue::YIon
			*/
			void getProtonDistribution(Map<UInt, double>& bb_charges, Map<UInt, double>& sc_charges, const AASequence& peptide, int charge,	Residue::ResidueType res_type = Residue::YIon);

			/** @brief calculates the charge state intensities of different charge states of the same ion
					
					@param peptide the peptide
					@param n_term_ion the prefix ion sequence
					@param c_term_ion the suffix ion sequence
					@param charge the charge
					@param n_term_type the ion type of the N-terminal ion; valid values are Residue::AIon, Residue::BIon
					@param n_term1 the probability of seeing a singly charged prefix ion 
					@param c_term1 the probability of seeing a singly charged suffix ion
					@param n_term2 the probability of seeing a doubly charged prefix ion
					@param c_term2 the probability of seeing a doubly charged suffix ion
					@param type the type of fragmentation (charge-directed, charge-remote of side chain)
			*/
			void getChargeStateIntensities(const AASequence& peptide, const AASequence& n_term_ion, const AASequence& c_term_ion, int charge, Residue::ResidueType n_term_type, double& n_term1,  double& c_term1, double& n_term2, double& c_term2, FragmentationType type);

			/// sets the proton distributions of the whole peptide, they are needed for the getChargeStateIntensities_ method and need to be recalculated each time if not given
			void setPeptideProtonDistribution(const Map<UInt, double>& bb_charge, const Map<UInt, double>& sc_charge);

			protected:

			// calculates the proton distribtion
			void calculateProtonDistribution_(const AASequence& peptide, int charge, Residue::ResidueType res_type = Residue::YIon, bool fixed_proton = false, UInt cleavage_site = 0, bool use_most_basic_site = false);
	
			// returns the proton affinity of the peptide with the given charge and ion type
			double getProtonAffinity_(const AASequence& ion, int charge, Residue::ResidueType res_type);

			// returns the (relative) Intensities of the possible charge states of the ion from peptide
			std::vector<double> getChargeStateIntensities_(const AASequence& peptide, const AASequence& ion, int charge, Residue::ResidueType res_type);

			// calculates the intensities of the different possible charge states
			void calcChargeStateIntensities_(const AASequence& peptide, const AASequence& n_term_ion, const AASequence& c_term_ion, int charge, Residue::ResidueType n_term_type,	double& n_term1, double& c_term1, double& n_term2, double& c_term2,	FragmentationType type);

			// returns the left and right GB values, NH2 and COOH if at terminus
			void getLeftAndRightGBValues_(const AASequence& peptide, double& left_gb, double& right_gb, UInt position);
			
			Map<UInt, double> sc_charge_;
			Map<UInt, double> bb_charge_;
			Map<UInt, double> sc_charge_full_;
			Map<UInt, double> bb_charge_full_;
			double E_;
			double E_c_term_;
			double E_n_term_;

	};
}
#endif // OPENMS_ANALYSIS_ID_PROTONDISTRIBUTIONMODEL_H

