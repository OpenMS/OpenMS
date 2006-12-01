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


#ifndef OPENMS_ANALYSIS_ID_PROTONDISTRIBUTIONMODEL_H
#define OPENMS_ANALYSIS_ID_PROTONDISTRIBUTIONMODEL_H

#include <OpenMS/DATASTRUCTURES/HashMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <vector>

namespace OpenMS 
{
	class ProtonDistributionModel 
	{
		public:
						
			ProtonDistributionModel();

			ProtonDistributionModel(const ProtonDistributionModel& model);
			
			virtual ~ProtonDistributionModel();

			ProtonDistributionModel& operator = (const ProtonDistributionModel& pdm);

      enum FragmentationType
      {
        ChargeDirected = 0,
        ChargeRemote,
        SideChain
      };

			void getProtonDistribution(	HashMap<Size, double>& bb_charges, 
																	HashMap<Size, double>& sc_charges, 
																	const AASequence& peptide,
																	int charge,
																	Residue::ResidueType res_type = Residue::YIon);

			void getChargeStateIntensities(	const AASequence& peptide,
																			const AASequence& n_term_ion,
																			const AASequence& c_term_ion,
																			int charge,
																			Residue::ResidueType n_term_type,
																			double& n_term1,
																			double& c_term1,
																			double& n_term2,
																			double& c_term2,
																			FragmentationType type);

			void setPeptideProtonDistribution(const HashMap<Size, double>& bb_charge, const HashMap<Size, double>& sc_charge);

			protected:

			void calculateProtonDistribution_(const AASequence& peptide, 
																				int charge, 
																				Residue::ResidueType res_type = Residue::YIon, 
																				bool fixed_proton = false, 
																				Size cleavage_site = 0,
																				bool use_most_basic_site = false);
	
			double getProtonAffinity_(const AASequence& ion, int charge, Residue::ResidueType res_type);

			// returns the (relative) Intensities of the possible charge states of the ion from peptide
			std::vector<double> getChargeStateIntensities_(const AASequence& peptide, const AASequence& ion, int charge, Residue::ResidueType res_type);

			void calcChargeStateIntensities_( const AASequence& peptide, 
																				const AASequence& n_term_ion,
																				const AASequence& c_term_ion,
																				int charge, 
																				Residue::ResidueType n_term_type,
																				double& n_term1,
																				double& c_term1,
																				double& n_term2,
																				double& c_term2,
																				FragmentationType type);

			void init_();

			HashMap<Size, double> sc_charge_;
			HashMap<Size, double> bb_charge_;
			HashMap<Size, double> sc_charge_full_;
			HashMap<Size, double> bb_charge_full_;

			HashMap<String, double> gb_sc_;
			HashMap<String, double> gb_bb_l_;
			HashMap<String, double> gb_bb_r_;

			double E_;
			double E_c_term_;
			double E_n_term_;

	};
}
#endif // OPENMS_ANALYSIS_ID_PROTONDISTRIBUTIONMODEL_H

