// -*- mode: C++; tab-width: 2; -*-
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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <gsl/gsl_randist.h>


#include <cmath>
#include <numeric>
#include <cstdlib>

#define COULOMB_REPULSION (DoubleReal)47.0  // from zhang: 47.0 kJ/mol
#define COULOMB_REPULSION2 (DoubleReal)47.0 // new

//#define CALC_CHARGE_STATES_DEBUG

using namespace std;

namespace OpenMS 
{
	ProtonDistributionModel::ProtonDistributionModel()
		:	DefaultParamHandler("ProtonDistributionModel"),
			E_(0),
			E_c_term_(0),
			E_n_term_(0)
	{
		defaults_.setValue("gb_bb_l_NH2", 916.84, "Gas-phase basicity value of N-terminus", StringList::create("advanced"));
		defaults_.setValue("gb_bb_r_COOH", -95.82, "Gas-phase basicity value of C-terminus", StringList::create("advanced"));
		defaults_.setValue("gb_bb_r_b-ion", 36.46, "Gas-phase basicity value of b-ion C-terminus", StringList::create("advanced"));
		defaults_.setValue("gb_bb_r_a-ion", 46.85, "Gas-phase basicity value of a-ion C-terminus", StringList::create("advanced"));
		defaults_.setValue("sigma", 0.5, "Width of the gaussian which distributes the mobile protons over the charge states, only for z > 3.", StringList::create("advanced"));
		defaults_.setValue("temperature", 500.0, "Temperature term ", StringList::create("advanced"));

		defaultsToParam_();
	}

	ProtonDistributionModel::~ProtonDistributionModel()
	{
	}

	ProtonDistributionModel::ProtonDistributionModel(const ProtonDistributionModel& model)
		: DefaultParamHandler(model),
			sc_charge_(model.sc_charge_),
	    bb_charge_(model.bb_charge_),
			sc_charge_full_(model.sc_charge_full_),
			bb_charge_full_(model.bb_charge_full_),
			E_(model.E_),
			E_c_term_(model.E_c_term_),
			E_n_term_(model.E_n_term_)
	{
	}
	
	ProtonDistributionModel& ProtonDistributionModel::operator = (const ProtonDistributionModel& model)
	{
		if (this != &model)
		{
			DefaultParamHandler::operator = (model);
			sc_charge_ = model.sc_charge_;
      bb_charge_  = model.bb_charge_;
      sc_charge_full_ = model.sc_charge_full_;
      bb_charge_full_ = model.bb_charge_full_;
      E_ = model.E_;
      E_c_term_ = model.E_c_term_;
      E_n_term_ = model.E_n_term_;
		}
		return *this;
	}
	
	void ProtonDistributionModel::setPeptideProtonDistribution(const vector<DoubleReal>& bb_charge, const vector<DoubleReal>& sc_charge)
	{
		bb_charge_full_ = bb_charge;
		sc_charge_full_ = sc_charge;
	}

	void ProtonDistributionModel::getProtonDistribution(vector<DoubleReal>& bb_charges,
															vector<DoubleReal>& sc_charges,
															const AASequence& peptide,
															Int charge,
															Residue::ResidueType res_type)
	{
		bb_charge_ = vector<DoubleReal>(peptide.size() + 1, 0.0);
		sc_charge_ = vector<DoubleReal>(peptide.size(), 0.0);
		calculateProtonDistribution_(peptide, charge, res_type);
		bb_charges = bb_charge_;
		sc_charges = sc_charge_;
	}

/*
	void ProtonDistributionModel::getChargeStateIntensities(const AASequence& peptide,
													                                const AASequence& n_term_ion,
																													const AASequence& c_term_ion,
																													Int charge,
																													Residue::ResidueType n_term_type,
																													DoubleReal& n_term1,
																													DoubleReal& c_term1,
																													DoubleReal& n_term2,
																													DoubleReal& c_term2,
																													FragmentationType type)
	{
		
		//calcChargeStateIntensities_(peptide, n_term_ion, c_term_ion, charge, n_term_type, n_term1, c_term1, n_term2, c_term2, type);
		vector<DoubleReal> n_term_intensities, c_term_intensities;
		calcChargeStateIntensities_(peptide, n_term_ion, c_term_ion, charge, n_term_type, n_term_intensities, c_term_intensities, type);
		n_term1 = n_term_intensities[0];
		c_term1 = c_term_intensities[0];
		if (charge == 1)
		{
			n_term2 = 0;
			c_term2 = 0;
		}
		else
		{
			n_term2 = n_term_intensities[1];
			c_term2 = c_term_intensities[1];
		}
		return;
	}
*/

	// this method calculates the distribution of protons, by using two separat ions given by the precursor peptide ion and the cleavage site
	// this is needed for charge directed cases, where the activated proton needs to be distributed between the two ions
	// cleavage_site is one-based; so b1/yk-1 would be cleavage_site == 1; it can be in the range [0,peptide.size()]
	void ProtonDistributionModel::calculateProtonDistributionIonPair_(const AASequence& peptide,
																																		Residue::ResidueType n_res_type,
																																		Size cleavage_site)
	{
		// TODO model this using one calculation for both ions

		DoubleReal q(0); // Zustandsumme

	  DoubleReal gb_bb_l_NH2 = (DoubleReal)param_.getValue("gb_bb_l_NH2");
	  DoubleReal gb_bb_r_COOH = (DoubleReal)param_.getValue("gb_bb_r_COOH");
	  DoubleReal gb_bb_r_bion = (DoubleReal)param_.getValue("gb_bb_r_b-ion");
  	DoubleReal gb_bb_r_aion = (DoubleReal)param_.getValue("gb_bb_r_a-ion");
		DoubleReal T = (DoubleReal)param_.getValue("temperature");

		// we calculate the distribution of only the last proton, all other protons are already distributed

		// so, first calculate the Zustandssumme Q of the N-term ion
		for (Size i = 0; i != cleavage_site; ++i)
  	{
    	// backbone energy
    	if (i == 0)
    	{
      	DoubleReal E = -(gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight());
      	q += (1.0 - bb_charge_[i]) * exp(-E * 1000 / (Constants::R * T));
    	}
    	else
    	{
      	if (i == cleavage_site - 1)
      	{
        	// position at the C-terminal end of the ion
        	DoubleReal E(-peptide[i].getBackboneBasicityLeft());
        	if (n_res_type == Residue::BIon)
        	{
          	E -= gb_bb_r_bion;
        	}
        	else
        	{
						// a-ion
            E -= gb_bb_r_aion;
          }
        	q += (1.0 - bb_charge_[i + 1]) * exp(-E * 1000 / (Constants::R * T));
				}

        // normal internal backbone position
        DoubleReal E = -(peptide[i - 1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
        q += (1.0 - bb_charge_[i]) * exp(-E * 1000 / (Constants::R * T));
    	}

   		// side chains
    	if (peptide[i].getSideChainBasicity() != 0)
    	{
      	DoubleReal E = -peptide[i].getSideChainBasicity();
      	q += (1.0 - sc_charge_[i]) * exp(-E * 1000 / (Constants::R * T));
    	}
  	}

		//cerr << "Q-N-term=" << 	q << endl;

		// add the parts of the C-term ion to the Zustandssumme
    for (Size i = cleavage_site; i != peptide.size(); ++i)
    {
      // backbone energy
      if (i == cleavage_site)
      {
        DoubleReal E = -(gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight());
				//cerr << "H3N: " << (1.0 - bb_charge_[i]) * exp(-E * 1000 / (Constants::R * T)) << endl;
        q += (1.0 - bb_charge_[i]) * exp(-E * 1000 / (Constants::R * T));
      }
      else
      {
        if (i == peptide.size() - 1)
        {
          // position at the C-terminal end of the ion, C-term ion always has COOH ending...
          DoubleReal E(0);
          E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_COOH);
					//cerr << "COOH: " << (1.0 - bb_charge_[i + 1]) * exp(-E * 1000 / (Constants::R * T)) << endl;
          q += (1.0 - bb_charge_[i + 1]) * exp(-E * 1000 / (Constants::R * T));
				}

				// normal internal backbone position
        DoubleReal E = -(peptide[i - 1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
				//cerr << "Internal: " << (1.0 - bb_charge_[i]) * exp(-E * 1000/ (Constants::R * T)) << endl;
        q += (1.0 - bb_charge_[i]) * exp(-E * 1000/ (Constants::R * T));
      }
      // side chains
      if (peptide[i].getSideChainBasicity() != 0)
      {
        DoubleReal E = -peptide[i].getSideChainBasicity();
				//cerr << "SideChain: " << (1.0 - sc_charge_[i]) * exp(-E * 1000 / (Constants::R * T)) << endl;
        q += (1.0 - sc_charge_[i]) * exp(-E * 1000 / (Constants::R * T));
      }
    }		

		//cerr << "Q-C-term=" << q << endl;
		
		// calculate proton availabilities of the N-terminal ion
		for (Size i = 0; i != cleavage_site; ++i)
  	{
    	// backbone
    	if (i == 0)
    	{
      	DoubleReal E = -(gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight());
      	bb_charge_ion_n_term_[i] = (1.0 - bb_charge_[i]) * exp(-E * 1000 /(Constants::R * T))/q;
    	}
    	else
    	{
      	if (i == cleavage_site - 1)
      	{
        	DoubleReal E(-peptide[i].getBackboneBasicityLeft());

        	if (n_res_type == Residue::BIon)
        	{
          	E -= gb_bb_r_bion;
        	}
        	else
        	{
						// a-ion
            E -= gb_bb_r_aion;
        	}
        	bb_charge_ion_n_term_[i + 1] = (1.0 - bb_charge_[i + 1]) * exp(-E * 1000/(Constants::R * T))/q;
				}

        // normal backbone position
        DoubleReal E = -(peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
        bb_charge_ion_n_term_[i] = (1.0 - bb_charge_[i]) * exp(-E * 1000 /(Constants::R * T))/q;
    	}

    	// side chains
    	if (peptide[i].getSideChainBasicity() != 0)
    	{
      	DoubleReal E = -peptide[i].getSideChainBasicity();
      	sc_charge_ion_n_term_[i] = (1.0 - sc_charge_[i]) * exp(-E * 1000 / (Constants::R * T))/q;
    	}
		}

		// same for the C-term ion
    for (Size i = cleavage_site; i != peptide.size(); ++i)
    {
      // backbone
      if (i == cleavage_site)
      {
        DoubleReal E = -(gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight());
        bb_charge_ion_c_term_[i - cleavage_site] = (1.0 - bb_charge_[i]) * exp(-E * 1000 /(Constants::R * T))/q;
      }
      else
      {
        if (i == peptide.size() - 1)
        {
          DoubleReal E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_COOH);
          bb_charge_ion_c_term_[i + 1 - cleavage_site] = (1.0 - bb_charge_[i + 1]) * exp(-E * 1000/(Constants::R * T))/q;
        }

        // normal backbone position
        DoubleReal E = -(peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
        bb_charge_ion_c_term_[i - cleavage_site] = (1.0 - bb_charge_[i]) * exp(-E * 1000 /(Constants::R * T))/q;
				//cerr << "Int: " << peptide[i].getOneLetterCode() << " " << bb_charge_[i] << 
      }

      // side chains
      if (peptide[i].getSideChainBasicity() != 0)
      {
        DoubleReal E = -peptide[i].getSideChainBasicity();
        sc_charge_ion_c_term_[i - cleavage_site] = (1.0 - sc_charge_[i]) * exp(-E * 1000 / (Constants::R * T))/q;
				//cerr << "SC: " << peptide[i].getOneLetterCode() << " " << sc_charge_[i] << " " << sc_charge_ion_c_term_[i - cleavage_site] << endl;
      }
    }

		return;
	}
	
	void ProtonDistributionModel::calculateProtonDistributionGreater2_(const AASequence& peptide, Int charge, Residue::ResidueType /*res_type*/)
	{
		//cerr << "void ProtonDistributionModel::calculateProtonDistributionGreater2_(" << peptide << ", " << charge << ", res_type=" << res_type << ")" << endl;
		//DoubleReal gb_bb_l_NH2 = param_.getValue("gb_bb_l_NH2");
    //DoubleReal gb_bb_r_COOH = param_.getValue("gb_bb_r_COOH");
    //DoubleReal gb_bb_r_bion = param_.getValue("gb_bb_r_b-ion");
    //DoubleReal gb_bb_r_aion = param_.getValue("gb_bb_r_a-ion");
		//const DoubleReal T(500.0);
					
		vector<DoubleReal> gb_bb, gb_sc;
		DoubleReal gb_left_n_term(0), gb_right_n_term(0);
		getLeftAndRightGBValues_(peptide, gb_left_n_term, gb_right_n_term, 0);
		gb_bb.push_back(gb_left_n_term + gb_right_n_term);
		Size count(1);
		for (AASequence::ConstIterator it = peptide.begin(); it != peptide.end(); ++it, ++count)
		{
			DoubleReal gb_left(0), gb_right(0);
			getLeftAndRightGBValues_(peptide, gb_left, gb_right, count);
			
			DoubleReal gb = gb_left + gb_right;
			gb_bb.push_back(gb);

			gb_sc.push_back(it->getSideChainBasicity());
		}
		
		// now distribute the charges until no site has more than 1.0 proton
		vector<DoubleReal> bb_coulomb(peptide.size() + 1, 0.0), sc_coulomb(peptide.size(), 0.0);
		Int actual_charge(charge);
		set<Size> sc_sites, bb_sites;
		while (true)
		{
			//cerr << "#protons remaining: " << actual_charge << endl;
			vector<DoubleReal> k_bb(peptide.size() + 1, 0.0), k_sc(peptide.size(), 0.0);
			count = 0;
			DoubleReal sum_k(0);
			for (vector<DoubleReal>::const_iterator it = gb_bb.begin(); it != gb_bb.end(); ++it, ++count)
			{
				if (bb_sites.find(count) == bb_sites.end())
				{
					k_bb[count] = exp((*it - bb_coulomb[count]) * 1000.0 / Constants::R / 500.0);
					sum_k += k_bb[count];
					//cerr << k_bb[count] << endl;
				}
			}

			count = 0;
			for (vector<DoubleReal>::const_iterator it = gb_sc.begin(); it != gb_sc.end(); ++it, ++count)
			{
				if (sc_sites.find(count) == sc_sites.end())
				{
					k_sc[count] = exp((*it - sc_coulomb[count]) * 1000.0 / Constants::R / 500.0);
					sum_k += k_sc[count];
					//cerr << k_sc[count] << endl;
				}
			}

			//cerr << "sum_k: " << sum_k << endl;

			vector<DoubleReal> p_bb(peptide.size() + 1, 1.0), p_sc(peptide.size(), 1.0);
			count = 0;
			for (vector<DoubleReal>::const_iterator it = k_bb.begin(); it != k_bb.end(); ++it, ++count)
			{
				if (bb_sites.find(count) == bb_sites.end())
				{
					p_bb[count] = (DoubleReal)actual_charge * *it / sum_k;
					bb_charge_[count] = p_bb[count];
				}
				//cerr << "BB" << count << ": " << p_bb[count] << endl;
			}

			count = 0;
			for (vector<DoubleReal>::const_iterator it = k_sc.begin(); it != k_sc.end(); ++it, ++count)
			{
				if (sc_sites.find(count) == sc_sites.end())
				{
					p_sc[count] = (DoubleReal)actual_charge * *it / sum_k;
					sc_charge_[count] = p_sc[count];
				}
				//cerr << "SC" << count << ": " << p_sc[count] << endl;
			}
	
			// check if there is a site containing more than one proton
			for (Size i = 0; i != p_bb.size(); ++i)
			{
				if (p_bb[i] > 1.0)
				{
					bb_charge_[i] = 1.0;
					bb_sites.insert(i);
					//cerr << "BackbonePosition " << i << " has charge > 1" << endl;
					--actual_charge;
				}
			}
			for (Size i = 0; i != p_sc.size(); ++i)
			{
				if (p_sc[i] > 1.0)
				{
					sc_charge_[i] = 1.0;
					//cerr << "SideChain " << i << " has charge > 1 " << endl;
					sc_sites.insert(i);
					--actual_charge;
				}
			}

			// now calculate the coloumb repulsions
			for (Size i = 0; i != gb_bb.size(); ++i)
			{
				// check if the site is not occupied by a "complete" proton
				if (bb_sites.find(i) == bb_sites.end())
				{
					DoubleReal coulomb_sum(0);
					for (set<Size>::const_iterator it = bb_sites.begin(); it != bb_sites.end(); ++it)
					{
						// calculate the distance between occupied site and this backbone site
						Size pos = *it;
						Size diff = (pos > i) ? pos - i : i - pos;
						coulomb_sum += COULOMB_REPULSION2 / (DoubleReal)diff;
					}

					for (set<Size>::const_iterator it = sc_sites.begin(); it != sc_sites.end(); ++it)
					{
						// calculate the distance between occupied side chain and this backbone site
						Size pos = *it;
						Size diff = (pos > i) ? pos -i : i - pos;
						++diff; // bond to the side chain counts extra
						coulomb_sum += COULOMB_REPULSION2 / (DoubleReal)diff;
					}
					bb_coulomb[i] = coulomb_sum;
					//cerr << "BB coulomb" << i << ": " << coulomb_sum << endl;
				}
			}
			for (Size i = 0; i != gb_sc.size(); ++i)
			{
				if (sc_sites.find(i) == sc_sites.end())
				{
					DoubleReal coulomb_sum(0);
					for (set<Size>::const_iterator it = bb_sites.begin(); it != bb_sites.end(); ++it)
					{
						Size pos = *it;
						Size diff = (pos > i) ? pos - i : i - pos;
						++diff;
						coulomb_sum += COULOMB_REPULSION2 / (DoubleReal)diff;
					}
					for (set<Size>::const_iterator it = sc_sites.begin(); it != sc_sites.end(); ++it)
					{
						Size pos = *it;
						Size diff = (pos > i) ? pos - i : i - pos;
						diff += 2;
						coulomb_sum += COULOMB_REPULSION2 / (DoubleReal)diff;
					}
					//cerr << "SC coulomb" << i << ": " << coulomb_sum << endl; 
					sc_coulomb[i] = coulomb_sum;
				}
			}
		
			// TODO think about what happens if #protons are greater than number of sites?!?
			if (bb_sites.size() == 0 && sc_sites.size() == 0)
			{
				break;
			}

			// search for entries > 1
			bool has_greater_one(false);
			for (vector<DoubleReal>::const_iterator it = p_bb.begin(); it != p_bb.end(); ++it)
			{
				if (*it > 1.0)
				{
					has_greater_one = true;
				}
			}

			for (vector<DoubleReal>::const_iterator it = p_sc.begin(); it != p_sc.end(); ++it)
			{
				if (*it > 1.0)
				{
					has_greater_one = true;
				}
			}

			if (!has_greater_one)
			{
				//cerr << "Has no site with more than 1.0 proton" << endl;
				break;
			}
		}
		
		return;

	}

	void ProtonDistributionModel::calculateProtonDistributionCharge2_(const AASequence& peptide,
																																		Residue::ResidueType res_type,
																																		bool fixed_proton,
																																		Size cleavage_site,
																																		bool use_most_basic_site)
	{
		DoubleReal q(0), sum_E(0), sum_E_n_term(0), sum_E_c_term(0); // Zustandsumme
		Size most_basic_site(0);
		bool most_basic_site_sc(false);

    DoubleReal gb_bb_l_NH2 = (DoubleReal)param_.getValue("gb_bb_l_NH2");
		DoubleReal gb_bb_r_COOH = (DoubleReal)param_.getValue("gb_bb_r_COOH");
		DoubleReal gb_bb_r_bion = (DoubleReal)param_.getValue("gb_bb_r_b-ion");
		DoubleReal gb_bb_r_aion = (DoubleReal)param_.getValue("gb_bb_r_a-ion");
		DoubleReal T = (DoubleReal)param_.getValue("temperature");

		if (!use_most_basic_site)
		{
			bb_charge_ = vector<DoubleReal>(peptide.size() + 1, 0.0);
			sc_charge_ = vector<DoubleReal>(peptide.size(), 0.0);
		}
		else
		{
			// find the most basic site
			DoubleReal max_prob(0);
			//cerr << "bb: ";
			for (Size i = 0; i != bb_charge_.size(); ++i)
			{
				//cerr << i << ". " << bb_charge_[i] << "; " << endl;
				if (bb_charge_[i] > max_prob)
				{
					max_prob = bb_charge_[i];
					most_basic_site = i;
				}
			}

			//cerr << endl << "sc: ";
			for (Size i = 0; i != sc_charge_.size(); ++i)
			{
				//cerr << i << ". " << sc_charge_[i] << "; " << endl;
				if (sc_charge_[i] > max_prob)
				{
					max_prob = sc_charge_[i];
					most_basic_site = i;
					most_basic_site_sc = true;
				}
			}
			//cerr << endl;
			

			bb_charge_ = vector<DoubleReal>(peptide.size() + 1, 0.0);
			sc_charge_ = vector<DoubleReal>(peptide.size(), 0.0);
		}

		Size fixed_site(0);
		if (fixed_proton)
		{
			fixed_site = cleavage_site;
		}

		bool fixed_site_sc(false);
		if (use_most_basic_site)
		{
			fixed_site = most_basic_site;
			fixed_site_sc = most_basic_site_sc;
		}

		for (Size i = 0; i != sc_charge_.size(); ++i)
		{
			sc_charge_[i] = 0;
		}

		for (Size i = 0; i != bb_charge_.size(); ++i)
		{
			bb_charge_[i] = 0;
		}
		//bb_charge_[peptide.size()] = 0;


	// fixed proton
	//
	// if two protons are available one proton is kept at the cleavage site
	// this is needed for the N/C-terminal charge distribution calculation
	//
	// use the fixed proton with precalculated charges of the other proton
	// -> get distribtion of the prior fixed one

	// make only sense with charge = 2
	if (fixed_proton || use_most_basic_site)
	{
		//cerr << "fixed site: " << fixed_site << " " << fixed_site_sc << endl;
		// fixed proton at fixed_site
		q = 0;
		
		DoubleReal gb_j(0);
		if (!fixed_site_sc)
		{
			if (fixed_site == 0)
			{
				gb_j = gb_bb_l_NH2 + peptide[fixed_site].getBackboneBasicityRight();
			}
			else
			{
				gb_j  = peptide[fixed_site - 1].getBackboneBasicityLeft() + peptide[fixed_site].getBackboneBasicityRight();
			}
		}
		else
		{
			gb_j = peptide[fixed_site].getSideChainBasicity();
		}
		
		for (Size i = 0; i <= peptide.size(); ++i)
		{
			DoubleReal gb_i(0);
      // proton 1 at N-terminus
      if (i == 0 || (i == cleavage_site && use_most_basic_site))
      {
				gb_i = gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight();
      }
      else
      {
        // proton 1 at N-terminus
        if (i == peptide.size())
        {
					if (res_type == Residue::BIon)
					{
						gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_bion;
					}
					else
					{
						if (res_type == Residue::AIon)
						{
							gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_aion;
						}
						else
						{
							gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_COOH;
						}
					}
        }
        else
        {
          // proton 1 at backbone
					gb_i = peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight();
        }
      }

			if (!fixed_site_sc)
			{
				if (i != fixed_site)
				{
					Int r_ij(abs((Int)i - (Int)(fixed_site)));
					q += exp(-(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000 / (Constants::R * T) - 500);
	
					DoubleReal gb_i_sc(0);
					if (i != peptide.size())
					{
						if (peptide[i].getSideChainBasicity() != 0)
						{
							gb_i_sc = peptide[i].getSideChainBasicity();
							q += exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500);
						}
					}
				}
				else
				{
					// last chance: the proton i is located at side chain of cleavage site
					if (i != peptide.size())
					{
						if (peptide[i].getSideChainBasicity() != 0)
						{
							DoubleReal gb_i_sc = peptide[i].getSideChainBasicity();
							q += exp(-(-gb_i - gb_i_sc + COULOMB_REPULSION) * 1000 / (Constants::R * T) - 500);
						}
					}
				}
			}
			else
			{
				// fixed site at side chain
				//
				// first, one proton at BB
				Int r_ij = abs((Int)i - (Int)fixed_site);
				q += exp(-(-gb_i - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 / (Constants::R * T) - 500);

				// only side chain site different from fixed one
				if (i != fixed_site && i != peptide.size())
				{
					DoubleReal gb_i_sc(0);
					gb_i_sc = peptide[i].getSideChainBasicity();
					q += exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 2)) * 1000 /(Constants::R * T) - 500);
				}
			}
		}

		// calculate availablities
		for (Size i = 0; i <= peptide.size(); ++i)
		{
			DoubleReal gb_i(0);
			if (i == 0 || (i == cleavage_site && use_most_basic_site))
			{
				gb_i = gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight();
			}
			else
			{
				if (i == peptide.size())
				{
					if (res_type == Residue::BIon)
					{
						gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_bion;
					}
					else
					{
						if (res_type == Residue::AIon)
						{
							gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_aion;
						}
						else
						{
							gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_COOH;
						}
					}
				}
				else
				{
					gb_i = peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight();
				}
			}
			
			if (!fixed_site_sc)
			{
				if (i != fixed_site)
				{
					Int r_ij(abs((Int)i - (Int)(fixed_site)));
					DoubleReal prob = exp(-(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000 / (Constants::R * T) - 500)/q;
					bb_charge_[i] += prob;
	
					DoubleReal add_E = exp(gb_i * 1000 / Constants::R / T);
					if (i < fixed_site - 1)
					{
						sum_E_n_term += add_E;
					}
					else
					{
						sum_E_c_term += add_E;
					}

					DoubleReal gb_i_sc(0);
					if (i != peptide.size())
					{
						if (peptide[i].getSideChainBasicity() != 0)
						{
							gb_i_sc = peptide[i].getSideChainBasicity();
							DoubleReal prob = exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500)/q;
							sc_charge_[i] += prob;
	
							DoubleReal add_E = exp(gb_i_sc * 1000 / Constants::R / T);
							if (i < fixed_site - 1)
 		       		{
 	  	       		sum_E_n_term += add_E;
      	  		}
        			else
        			{
        	 	 		sum_E_c_term += add_E;
        			}
						}
					}	
				}
				else
				{
					if (i != peptide.size())
					{
						// SC position
						DoubleReal gb_i_sc(0);
						if (peptide[i].getSideChainBasicity() != 0)
						{
							gb_i_sc = peptide[i].getSideChainBasicity();
							DoubleReal prob = exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION) * 1000 /(Constants::R * T) - 500)/q;
							sc_charge_[i] += prob;

							DoubleReal add_E = exp(gb_i_sc * 1000 / Constants::R / T);
        			if (i < fixed_site - 1)
        			{
          			sum_E_n_term += add_E;
        			}
        			else
        			{
          			sum_E_c_term += add_E;
        			}
						}
					}
				}
			}
			else
			{
				// fixed site at side chain
				Int r_ij = abs((Int)i - (Int)fixed_site);
				DoubleReal prob = exp(-(-gb_i - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 / (Constants::R * T) - 500)/q;
				bb_charge_[i] += prob;

				DoubleReal add_E = exp(gb_i * 1000 / Constants::R / T);
				if (i <= fixed_site)
				{
					sum_E_n_term += add_E;
				}
				else
				{
					sum_E_c_term += add_E;
				}

				if (i != fixed_site && i != peptide.size())
				{
					DoubleReal gb_i_sc(0);
					if (peptide[i].getSideChainBasicity() != 0)
					{
						gb_i_sc = peptide[i].getSideChainBasicity();
						DoubleReal prob = exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 2)) * 1000 / (Constants::R * T) - 500)/q;
						sc_charge_[i] += prob;

						DoubleReal add_E = exp(gb_i_sc * 1000 / Constants::R / T);
						if (i <= fixed_site)
						{
							sum_E_n_term += add_E;
						}
						else
						{
							sum_E_c_term += add_E;
						}
					}
				}
			}
		}
	}

	
	// DoubleReal charged
	if (!fixed_proton && !use_most_basic_site)
	{
		// calculate sum
		Int count(0);
		for (Size i = 0; i <= peptide.size(); ++i)
		{
			for (Size j = i; j <= peptide.size(); ++j)
			{
				DoubleReal gb_i(0), gb_j(0);
				// proton 1 at N-terminus
				if (i == 0)
				{
					//gb_i = gb_bb_l_["NH2"] + gb_bb_r_[peptide[i].getOneLetterCode()];
					gb_i = gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight();
				}
				else
				{
					String aa_i_l(peptide[i-1].getOneLetterCode());
					// proton 1 at N-terminus
					if (i == peptide.size())
					{
						if (res_type == Residue::BIon)
						{
							//gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["b-ion"];
							gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_bion;
						}
						else
						{
							if (res_type == Residue::AIon)
							{
								//gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["a-ion"];
								gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_aion;
							}
							else
							{
								//gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["COOH"];
								gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_COOH;
							}
						}
					}
					else
					{
						// proton 1 at backbone
						//gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_[peptide[i].getOneLetterCode()];
						gb_i = peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight();
					}
				}
				// proton 2 at N-terminus
				if (j == 0)
				{
					//gb_j = gb_bb_l_["NH2"] + gb_bb_r_[peptide[j].getOneLetterCode()];
					gb_j = gb_bb_l_NH2 + peptide[j].getBackboneBasicityRight();
				}
				else
				{
					//String aa_j_l(peptide[j-1].getOneLetterCode());
					DoubleReal gb_j_l = peptide[j-1].getBackboneBasicityLeft();
					// proton 2 at C-terminus
					if (j == peptide.size())
					{
						if (res_type == Residue::BIon)
						{
							//gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["b-ion"];
							gb_j = gb_j_l + gb_bb_r_bion;
						}
						else
						{
							if (res_type == Residue::AIon)
							{
								//gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["a-ion"];
								gb_j = gb_j_l + gb_bb_r_aion;
							}
							else
							{
								//gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["COOH"];
								gb_j = gb_j_l + gb_bb_r_COOH;
							}
						}
					}
					else
					{
						// proton 2 at backbone
						//gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_[peptide[j].getOneLetterCode()];
						gb_j = gb_j_l + peptide[j].getBackboneBasicityRight();
					}
				}
				if (i != j)
				{
					// distance of protons
					Int r_ij(abs((Int)i - (Int)j));
					q += exp(-(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000 /(Constants::R * T) - 500);
					//cerr << "1.\t" << -(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000/(Constants::R * T) << endl;
					++count;

					DoubleReal gb_i_sc(0), gb_j_sc(0);
					if (i != peptide.size())
					{
						// side chain of proton 1
						if (/*gb_sc_.has(peptide[i].getOneLetterCode())*/peptide[i].getSideChainBasicity() != 0)
						{
							//gb_i_sc = gb_sc_[peptide[i].getOneLetterCode()];
							gb_i_sc = peptide[i].getSideChainBasicity();
							q += exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500);
							//cerr << "2.\t" << -(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000/(Constants::R * T) << endl;
							++count;
						}
					}
					if (j != peptide.size())
					{
						// side chain of proton 2
						if (/*gb_sc_.has(peptide[j].getOneLetterCode())*/peptide[j].getSideChainBasicity() != 0)
						{
							//gb_j_sc = gb_sc_[peptide[j].getOneLetterCode()];
							gb_j_sc = peptide[j].getSideChainBasicity();
							q += exp(-(-gb_i - gb_j_sc + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500);
							//cerr << "3.\t" << -(-gb_i - gb_j_sc + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500 << endl;
							++count;
							// both at side chain?
							if (gb_i_sc != 0)
							{
								q += exp(-(-gb_i_sc - gb_j_sc + COULOMB_REPULSION/(r_ij + 2)) * 1000/(Constants::R * T) - 500);
								//cerr << "4.\t" << -(-gb_i_sc - gb_j_sc + COULOMB_REPULSION/(r_ij + 2)) * 1000 /(Constants::R * T) - 500 << endl;
								++count;
							}
						}
					}
				}
				else
				{
					if (i != peptide.size())
					{
						// one at side chain, the other one at backbone of same residue
						if (/*gb_sc_.has(peptide[i].getOneLetterCode())*/peptide[i].getSideChainBasicity() != 0)
						{
							//DoubleReal gb_i_sc = gb_sc_[peptide[i].getOneLetterCode()];
							DoubleReal gb_i_sc = peptide[i].getSideChainBasicity();
							q += exp(-(-gb_i - gb_i_sc + COULOMB_REPULSION) * 1000 / (Constants::R * T) - 500);
							//cerr << "5.\t" << -(-gb_i - gb_i_sc + COULOMB_REPULSION) * 1000/ (Constants::R * T) -500 << endl;
							++count;
						}
					}
				}
			}
		}

		// calculate availabilities
		for (Size i = 0; i <= peptide.size(); ++i)
		{
			for (Size j = i; j <= peptide.size(); ++j)
			{
				DoubleReal gb_i(0), gb_j(0);
				// calculate the backbone proton gb's
				// N-terminus
				if (i == 0)
				{
					//gb_i = gb_bb_l_["NH2"] + gb_bb_r_[peptide[i].getOneLetterCode()];
					gb_i = gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight();
				}
				else
				{
					//String aa_i_l(peptide[i-1].getOneLetterCode());
					
					// C-terminus
					if (i == peptide.size())
					{
						if (res_type == Residue::BIon)
						{
							//gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["b-ion"];
							gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_bion;
						}
						else
						{
							if (res_type == Residue::AIon)
							{
								//gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["a-ion"];
								gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_aion;
							}
							else
							{
								//gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["COOH"];
								gb_i = peptide[i-1].getBackboneBasicityLeft() + gb_bb_r_COOH;
							}
						}
					}
					else
					{
						// internal BB gb's
						//gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_[peptide[i].getOneLetterCode()];
						gb_i = peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight();
					}
				}
				// N-terminus
				if (j == 0)
				{
					//gb_j = gb_bb_l_["NH2"] + gb_bb_r_[peptide[j].getOneLetterCode()];
					gb_j = gb_bb_l_NH2 + peptide[j].getBackboneBasicityRight();
				}
				else
				{
					String aa_j_l(peptide[j-1].getOneLetterCode());
					// C-terminus
					if (j == peptide.size())
					{
						if (res_type == Residue::BIon)
						{
							//gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["b-ion"];
							gb_j = peptide[j-1].getBackboneBasicityLeft() + gb_bb_r_bion;
						}
						else
						{
							if (res_type == Residue::AIon)
							{
								//gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["a-ion"];
								gb_j = peptide[j-1].getBackboneBasicityLeft() + gb_bb_r_aion;
							}
							else
							{
								//gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["COOH"];
								gb_j = peptide[j-1].getBackboneBasicityLeft() + gb_bb_r_COOH;
							}
						}
					}
					else
					{
						//gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_[peptide[j].getOneLetterCode()];
						gb_j = peptide[j-1].getBackboneBasicityLeft() + peptide[j].getBackboneBasicityRight();
					}
				}

				// protons at different residues
				if (i != j)
				{
					// distance of the protons
					Int r_ij(abs((Int)i - (Int)j));
					// calc probability
					DoubleReal prob = exp(-(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000 / (Constants::R * T) - 500)/q;
					// add prob to site of first proton
					bb_charge_[i] += prob;
					// add to apperent GB
					bb_charge_[j] += prob;

					// side chains
					DoubleReal gb_i_sc(0), gb_j_sc(0);
					if (i != peptide.size())
					{
						if (/*gb_sc_.has(peptide[i].getOneLetterCode())*/peptide[i].getSideChainBasicity() != 0)
						{
							//gb_i_sc = gb_sc_[peptide[i].getOneLetterCode()];
							gb_i_sc = peptide[i].getSideChainBasicity();
							DoubleReal prob = exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500)/q;
							sc_charge_[i] += prob;
							bb_charge_[j] += prob;
						}
					}
					
					if (j != peptide.size())
					{
						if (/*gb_sc_.has(peptide[j].getOneLetterCode())*/peptide[j].getSideChainBasicity() != 0)
						{
							//gb_j_sc = gb_sc_[peptide[j].getOneLetterCode()];
							gb_j_sc = peptide[j].getSideChainBasicity();
							DoubleReal prob = exp(-(-gb_i - gb_j_sc + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500)/q;
							bb_charge_[i] += prob;
							sc_charge_[j] += prob;

							// both protons at sidechains
							if (gb_i_sc != 0)
							{
								DoubleReal prob = exp(-(-gb_i_sc - gb_j_sc + COULOMB_REPULSION/(r_ij + 2)) * 1000 /(Constants::R * T) - 500)/q;
								sc_charge_[i] += prob;
								sc_charge_[j] += prob;
							}
						}
					}
				}
				else
				{
					// protons at the same residue
					if (i != peptide.size())
					{
						if (/*gb_sc_.has(peptide[i].getOneLetterCode())*/peptide[i].getSideChainBasicity() != 0)
						{
							//DoubleReal gb_i_sc = gb_sc_[peptide[i].getOneLetterCode()];
							DoubleReal gb_i_sc = peptide[i].getSideChainBasicity();
							DoubleReal prob = exp(-(-gb_i - gb_i_sc + COULOMB_REPULSION) * 1000 / (Constants::R * T) - 500)/q;
							sc_charge_[i] += prob;
							sc_charge_[j] += prob;
						}
					}
				}
			}
		}
	}	

	#if 0
	cerr << "side chain proton availabilities" << endl;
	DoubleReal sum(0);
	for (unsigned Int i = 0; i != peptide.size(); ++i)
	{
		if (sc_charge_.has(i))
		{
			cerr << i << ".\t" << peptide[i].getThreeLetterCode() << ": " << sc_charge_[i] << endl;
			sum += sc_charge_[i];
		}
		else
		{
			cerr << i << ".\t" << peptide[i].getThreeLetterCode() << ": 0" << endl;
		}
	}

	cerr << "\nbackbone proton availabilities" << endl;
	for (unsigned Int i = 0; i  <= peptide.size(); ++i)
	{
		if (bb_charge_.has(i))
		{
			cerr << i << ".\t" << bb_charge_[i] << endl;
			sum += bb_charge_[i];
		}
		else
		{
			cerr << i << ".\t" << 0 << endl;
		}
	}
	
	cerr << "(sum=" << sum << ")" << endl;
	#endif
	
	E_ = sum_E;
	if (fixed_proton)
	{
		E_n_term_ = sum_E_n_term;
		E_c_term_ = sum_E_c_term;
	}
	else
	{
		E_n_term_ = 0;
		E_c_term_ = 0;
	}

	}

	void ProtonDistributionModel::calculateProtonDistributionCharge1_(const AASequence& peptide, Residue::ResidueType res_type)
  {
	// single charged
	DoubleReal q(0), sum_E(0)/*, sum_E_n_term(0), sum_E_c_term(0)*/; // Zustandsumme

  DoubleReal gb_bb_l_NH2 = (DoubleReal)param_.getValue("gb_bb_l_NH2");
  DoubleReal gb_bb_r_COOH = (DoubleReal)param_.getValue("gb_bb_r_COOH");
  DoubleReal gb_bb_r_bion = (DoubleReal)param_.getValue("gb_bb_r_b-ion");
  DoubleReal gb_bb_r_aion = (DoubleReal)param_.getValue("gb_bb_r_a-ion");
	DoubleReal T = (DoubleReal)param_.getValue("temperature");

	for (Size i = 0; i != peptide.size(); ++i)
	{
		// backbone energy
		if (i == 0)
		{
			DoubleReal E = -(gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight());
			q += exp(-E * 1000 / (Constants::R * T));
		}
		else
		{
			if (i == peptide.size() - 1)
			{
				// position at the C-terminal end of the ion
				DoubleReal E(0);
				if (res_type == Residue::BIon)
				{
					E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_bion);
				}
				else
				{
					if (res_type == Residue::AIon)
					{
						E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_aion);
					}
					else
					{
						E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_COOH);
					}
				}
				q += exp(-E * 1000 / (Constants::R * T));
				E = -(peptide[i - 1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
				q += exp(-E * 1000/ (Constants::R * T));
			}
			else
			{
				// normal internal backbone position
				DoubleReal E = -(peptide[i - 1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
				q += exp(-E * 1000 / (Constants::R * T));
			}
		}
		
		// side chains
		if (peptide[i].getSideChainBasicity() != 0)
		{
			DoubleReal E = -peptide[i].getSideChainBasicity();
			q += exp(-E * 1000 / (Constants::R * T));
		}
	}

	#if 0
	cout << "Q=" << q << endl;
	#endif
	// calculate the availabilities
	for (Size i = 0; i != peptide.size(); ++i)
	{
		// backbone
		if (i == 0)
		{
			DoubleReal E = -(gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight());
			bb_charge_[i] = exp(-E * 1000 /(Constants::R * T))/q;
			sum_E += exp(-E * 1000/Constants::R/T);
		}
		else
		{
			if (i == peptide.size() - 1)
			{
				DoubleReal E(0);
				
				if (res_type == Residue::BIon)
				{
					E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_bion);
				}
				else
				{
					if (res_type == Residue::AIon)
					{
						E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_aion);
					}
					else
					{
						E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_COOH);
					}
				}
				// TODO charge order; bug????
				bb_charge_[i + 1] = exp(-E * 1000/(Constants::R * T))/q;
				sum_E += exp(-E * 1000/Constants::R/T);

				E = -(peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
				bb_charge_[i] = exp(-E * 1000 /(Constants::R * T))/q;
				sum_E += exp(-E * 1000/Constants::R/T);
			}
			else
			{
				// normal backbone position
				DoubleReal E = -(peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
				bb_charge_[i] = exp(-E * 1000 /(Constants::R * T))/q;
				sum_E += exp(-E * 1000/Constants::R/T);
			}
		}

		// side chains
		if (peptide[i].getSideChainBasicity() != 0)
		{
			DoubleReal E = -peptide[i].getSideChainBasicity();
			sc_charge_[i] = exp(-E * 1000 / (Constants::R * T))/q;
			sum_E += exp(-E * 1000/Constants::R/T);
		}
  }
		E_ = sum_E;
	}
	
	void ProtonDistributionModel::calculateProtonDistribution_(const AASequence& peptide, 
																								Int charge, Residue::ResidueType res_type, 
																								bool fixed_proton, 
																								Size cleavage_site,
																								bool use_most_basic_site)
	{
		if (charge == 1)
		{
			calculateProtonDistributionCharge1_(peptide, res_type);
			return;
		}
		if (charge == 2)
		{
			calculateProtonDistributionCharge2_(peptide, res_type, fixed_proton, cleavage_site, use_most_basic_site);
			return;
		}
		
		// charge > 2
		calculateProtonDistributionGreater2_(peptide, charge, res_type);
		return;
	}

/*
	DoubleReal ProtonDistributionModel::getProtonAffinity_(const AASequence& peptide, Int charge, Residue::ResidueType res_type)
	{
		//const DoubleReal T(500.0);

		DoubleReal pa(0);
		calculateProtonDistribution_(peptide, charge, res_type);

		//pa = Constants::R * T * log(E_);

		// new test
		DoubleReal sum(0);
		if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
		{
			for (Size i = 1; i <= peptide.size(); ++i)
			{
				sum += bb_charge_[i];
				if (sc_charge_.has(i-1))
				{
					sum += sc_charge_[i-1];
				}
			}
			sum += bb_charge_[0]/2;
		}
		else
		{
			if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
			{
				for (Size i = bb_charge_.size() - peptide.size() - 1; i != bb_charge_.size(); ++i)
				{
					sum += bb_charge_[i];
					if (sc_charge_.has(i))
					{
						sum += sc_charge_[i];
					}
				}
				sum += bb_charge_[0]/2;
			}
		}

		pa = sum;
	
		//return pa;
		return E_;
	}
*/

	void ProtonDistributionModel::getChargeStateIntensities(const AASequence& peptide, 
																 													const AASequence& n_term_ion, 
																													const AASequence& c_term_ion, 
																													Int charge, 
																													Residue::ResidueType n_term_type,
																													std::vector<DoubleReal>& n_term_intensities, 
																													std::vector<DoubleReal>& c_term_intensities, 
																													FragmentationType type)
	{
		calcChargeStateIntensities_(peptide, n_term_ion, c_term_ion, charge, n_term_type, n_term_intensities, c_term_intensities, type);
		return;
	}

	void ProtonDistributionModel::calcChargeStateIntensities_(const AASequence& peptide, 
																														const AASequence& n_term_ion, 
																														const AASequence& c_term_ion, 
																														Int charge,
																														Residue::ResidueType n_term_type,
																														vector<DoubleReal>& n_term_intensities,
																														vector<DoubleReal>& c_term_intensities,
																														FragmentationType type)
	{
		// original method works well
		if (charge == 1)
		{
			DoubleReal c_term_int1(0), c_term_int2(0), n_term_int1(0), n_term_int2(0);
			n_term_intensities.clear();
			c_term_intensities.clear();
			calcChargeStateIntensities_(peptide, n_term_ion, c_term_ion, charge, n_term_type, n_term_int1, c_term_int1, n_term_int2, c_term_int2, type);
			n_term_intensities.push_back(n_term_int1);
			c_term_intensities.push_back(c_term_int1);
			return;
		}

		if (charge == 2)
		{
			DoubleReal c_term_int1(0), c_term_int2(0), n_term_int1(0), n_term_int2(0);
      n_term_intensities.clear();
      c_term_intensities.clear();
      calcChargeStateIntensities_(peptide, n_term_ion, c_term_ion, charge, n_term_type, n_term_int1, c_term_int1, n_term_int2, c_term_int2, type);
      n_term_intensities.push_back(n_term_int1);
			n_term_intensities.push_back(n_term_int2);
      c_term_intensities.push_back(c_term_int1);
			c_term_intensities.push_back(c_term_int2);
			return;
		}
	

		// charge > 2
		n_term_intensities = vector<DoubleReal>(charge, 0.0);
		c_term_intensities = vector<DoubleReal>(charge, 0.0);
		
		// calculate the number of in-active protons
		// this is simply the number of protons in case of charge-remote and 
		// side-chain induced cleavages (side-chain protons stay at the side chain)
		// for charge-directed cleavages there must be one proton which induces 
		// the cleavage, however, this can be distributed over several places
		Int num_active_protons = charge;
		if (type == ChargeDirected)
		{
			num_active_protons = charge - 1;
		}
		
		// calc proton distribution
		calculateProtonDistribution_(peptide, num_active_protons, Residue::Full);
    DoubleReal n_term_sum(0), c_term_sum(0);

		// sum up all protons located at the N-term/C-term part of the peptide
    for (Size i = 0; i != n_term_ion.size(); ++i)
    {
			//cerr << "N-term: i=" << i << ", " << peptide[i].getOneLetterCode() << ", BB-charge=" << bb_charge_[i] << ", SC-charge=" << sc_charge_[i] << endl;
      n_term_sum += bb_charge_[i];
      n_term_sum += sc_charge_[i];
    }
    for (Size i = n_term_ion.size(); i != peptide.size(); ++i)
    {
			//cerr << "C-term: i=" << i << ", " << peptide[i].getOneLetterCode() << ", BB-charge=" << bb_charge_[i] << ", SC-charge=" << sc_charge_[i] << endl;
      c_term_sum += bb_charge_[i + 1];
      c_term_sum += sc_charge_[i];
    }
    //c_term_sum += bb_charge_[peptide.size()];
		//cerr << n_term_ion << " " << c_term_ion << " " << n_term_sum << " " << c_term_sum << endl;
		//n_term_intensities[0] = n_term_sum;
		//c_term_intensities[0] = c_term_sum;
		// now we have the distributions of the protons
		
		//cerr << "Init: " << n_term_sum << " " << c_term_sum << endl;
		if (type == ChargeDirected)
		{
			// charge directed case
			// the proton which induces the cleavage is handled separately
			bb_charge_ion_n_term_ = vector<DoubleReal>(n_term_ion.size() + 1, 0.0);
			bb_charge_ion_c_term_ = vector<DoubleReal>(c_term_ion.size() + 1, 0.0);
			sc_charge_ion_n_term_ = vector<DoubleReal>(n_term_ion.size(), 0.0);
			sc_charge_ion_c_term_ = vector<DoubleReal>(c_term_ion.size(), 0.0);
			calculateProtonDistributionIonPair_(peptide, n_term_type, n_term_ion.size());
			//cerr << "NTerm: ";
			for (Size i = 0; i != n_term_ion.size(); ++i)
			{
				//cerr << i << ", bb=" << bb_charge_ion_n_term_[i] << ", sc=" << sc_charge_ion_n_term_[i] << " ";
				n_term_sum += bb_charge_ion_n_term_[i];
				n_term_sum += sc_charge_ion_n_term_[i];
			}
			//cerr << bb_charge_ion_n_term_[n_term_ion.size()] << endl;
			n_term_sum += bb_charge_ion_n_term_[n_term_ion.size()];

			//cerr << "CTerm: ";
			for (Size i = 0; i != c_term_ion.size(); ++i)
			{
				//cerr << i << ", bb=" << bb_charge_ion_c_term_[i] << ", sc=" << sc_charge_ion_c_term_[i] << " ";
				c_term_sum += bb_charge_ion_c_term_[i];
				c_term_sum += sc_charge_ion_c_term_[i];
			}
			//cerr << bb_charge_ion_c_term_[c_term_ion.size()] << endl;
			c_term_sum += bb_charge_ion_c_term_[c_term_ion.size()];
		}

    // we simply need to calculate the charge state distribution according
    // to the proton probabilities we calculated above
		DoubleReal sigma = (DoubleReal)param_.getValue("sigma");
    for (Int z = 1; z <= charge; ++z)
    {
      n_term_intensities[z - 1] = gsl_ran_gaussian_pdf(fabs(n_term_sum - (DoubleReal)z), sigma);
      c_term_intensities[z - 1] = gsl_ran_gaussian_pdf(fabs(c_term_sum - (DoubleReal)z), sigma);
    }

		return;
	}

	void ProtonDistributionModel::calcChargeStateIntensities_( const AASequence& peptide, 
																									const AASequence& n_term_ion, 
																									const AASequence& c_term_ion,
																									Int charge, 
																									Residue::ResidueType n_term_type,	
																									DoubleReal& n_term1, 
																									DoubleReal& c_term1, 
																									DoubleReal& n_term2, 
																									DoubleReal& c_term2,
																									FragmentationType type)
	{
	
		DoubleReal n_term_kapp(0), c_term_kapp(0);
		if (charge == 1)
		{
			if (type == ChargeDirected || type == ChargeRemote)
			{
				// get the K_app of N and C-terminal fragment respectively
				calculateProtonDistribution_(n_term_ion, 1, n_term_type);
				n_term_kapp = E_;
				calculateProtonDistribution_(c_term_ion, 1, Residue::YIon);
				c_term_kapp = E_;

				// calc the ratio
				n_term1 = n_term_kapp / (n_term_kapp + c_term_kapp);
				c_term1 = c_term_kapp / (n_term_kapp + c_term_kapp);
				//
				
				//DoubleReal pa_n = log(n_term_kapp);
				//DoubleReal pa_c = log(c_term_kapp);

				//DoubleReal ratio_bx_yz = exp(pa_n - pa_c);

				//n_term1 = ratio_bx_yz/*ratio_bx_yz / (1.0 + ratio_bx_yz)*/;
				//c_term1 = 1.0 /*/ (ratio_bx_yz + 1.0)*/;

				// of course ++ ions are not available
				n_term2 = 0;
				c_term2 = 0;
				

				//cerr << "ChargeStateIntensities: " << n_term_ion << " - " << c_term_ion << " z=1 " << n_term_kapp << " " << c_term_kapp << " " << n_term1 << " " << c_term1 << endl;
			}
			else
			{
				if (type == SideChain)
				{
					// the proton stays at the fragmentation site (N-terminal fragment)
					n_term1 = 1;
					c_term1 = 0;
					n_term2 = 0;
					c_term2 = 0;
				}
				else
				{
					// not possible
					cerr << "calcChargeStateIntensities_: unknown fragmentation type (" << type << ")" << endl;
				}
			}
			return;
		}

		if (charge == 2)
		{
			if (type == ChargeDirected)
			{
				// calc proton distribution with one fixed at cleavage site
				calculateProtonDistribution_(peptide, 2, Residue::Full, true, n_term_ion.size());
				//calculateProtonDistribution_(peptide, 1, Residue::Full);
				DoubleReal p_n(0), p_c(0);

				p_n = E_n_term_ / (E_n_term_ + E_c_term_);
				if (p_n < 0)
				{
					p_n = 0;
				}
				p_c = E_c_term_ / (E_n_term_ + E_c_term_);
				if (p_c < 0)
				{
					p_c = 0;
				}

#ifdef CALC_CHARGE_STATES_DEBUG
				cerr << "E_n_term_=" << E_n_term_ << ", E_c_term_=" << E_c_term_ << ", p_n=" << p_n << ", p_c=" << p_c << endl;
#endif
								
	
				// calc proton distribution of second proton with other one at most basic site fixed
				calculateProtonDistribution_(peptide, 2, Residue::Full, false, n_term_ion.size(), true);

#ifdef CALC_CHARGE_STATES_DEBUG
				cerr << "Distribution of second proton: " << endl;
				cerr << "BB: ";
				for (Size i = 0; i != bb_charge_.size(); ++i)
				{
					cerr << "; " << i << ". " << bb_charge_[i];
				}
				cerr << "\nSC: ";
				for (Size i = 0; i != sc_charge_.size(); ++i)
				{
					cerr << "; " << i << ". " << sc_charge_[i];
				}
				cerr << endl;
#endif
				
				DoubleReal singly_charged(0);
				for (Size i = 0; i != n_term_ion.size(); ++i)
				{
					n_term2 += bb_charge_[i] * p_n;
					singly_charged += bb_charge_[i] * p_c;
					if (sc_charge_[i] != 0)
					{
						n_term2 += sc_charge_[i] * p_n;
						singly_charged += sc_charge_[i]  * p_c;
					}
				}

				for (Size i = n_term_ion.size(); i <= peptide.size(); ++i)
				{
					c_term2 += bb_charge_[i] * p_c;
					singly_charged += bb_charge_[i] * p_n;
					if (i < peptide.size() && sc_charge_[i] != 0)
					{
						c_term2 += sc_charge_[i]  * p_c;
						singly_charged += sc_charge_[i] * p_n;
					}
				}

				n_term1 = singly_charged;
				c_term1 = singly_charged;

				//cerr << E_n_term_ << "\t" << E_c_term_ << "\t" << p_n << "\t" << p_c << "\t" << endl;


				// calculate charge losses
				DoubleReal gb_n_term = AAIndex::calculateGB(n_term_ion);
				DoubleReal gb_c_term = AAIndex::calculateGB(c_term_ion);

				DoubleReal b(828.18); // kJ/mol

				DoubleReal gb_n_term_loss = exp(- (gb_n_term - b)/1000.0);
				DoubleReal gb_c_term_loss = exp(- (gb_c_term - b)/1000.0);

#ifdef CALC_CHARGE_STATES_DEBUG
				cerr << "Loss: N-term: " << gb_n_term << " -> " << gb_n_term_loss << ", C-term: " << gb_c_term << " -> " << gb_c_term_loss << endl;
#endif

				n_term1 += n_term2 * (1.0 - gb_n_term_loss);
				n_term2 *= gb_n_term_loss;
				c_term1 += c_term2 * (1.0 - gb_c_term_loss);
				c_term2 *= gb_c_term_loss;


				// TODO normalization correct?
				DoubleReal sum(0);
				sum += n_term1 + n_term2 + c_term1 + c_term2;
				n_term1 /= sum;
				n_term2 /= sum;
				c_term1 /= sum;
				c_term2 /= sum;

#ifdef CALC_CHARGE_STATES_DEBUG
				cerr << "CD:     charge=2, " << n_term_ion << "|" << c_term_ion << ", n_term1=" << n_term1 << ", n_term2=" << n_term2 << ", c_term1=" << c_term1 << ", c_term2=" << c_term2 << endl;
#endif
			}
			else
			{
				if (type == ChargeRemote || type == SideChain)
				{
					// TODO ranges correct? Missing some sites of the peptide
					DoubleReal n_term_sum(0), c_term_sum(0);
					for (Size i = 0; i != n_term_ion.size(); ++i)
					{
						n_term_sum += bb_charge_full_[i];
						n_term_sum += sc_charge_full_[i];
					}
					for (Size i = n_term_ion.size(); i != peptide.size(); ++i)
					{
						c_term_sum += bb_charge_full_[i];
						c_term_sum += sc_charge_full_[i];
					}
					c_term_sum += bb_charge_full_[peptide.size()];
						
					if (n_term_sum - 1 > 0)
					{
						n_term2 = n_term_sum - 1;
						n_term1 = 1 - n_term2;
					}
					else
					{
						n_term1 = n_term_sum;
						n_term2 = 0;
					}
					if (c_term_sum - 1 > 0)
					{
						c_term2 = c_term_sum - 1;
						c_term1 = 1 - c_term2;
					}
					else
					{
						c_term1 = c_term_sum;
						c_term2 = 0;
					}

					DoubleReal sum(0);
					sum += n_term1 + n_term2 + c_term1 + c_term2;
					n_term1 /= sum;
					n_term2 /= sum;
					c_term1 /= sum;
					c_term2 /= sum;
#ifdef CALC_CHARGE_STATES_DEBUG
					cerr << "CR/SC: charge=2, " << n_term_ion << "|" << c_term_ion << ", n_term1=" << n_term1 << ", n_term2=" << n_term2 << ", c_term1=" << c_term1 << ", c_term2=" << c_term2 << endl;
#endif
				}
				else
				{
					cerr << "calcChargeStateIntensities_: unknown fragmentation type (" << type << ")" << endl;
				}
			}
		}
		if (charge > 2)
		{
			/*const AASequence& peptide, const AASequence& n_term_ion, const AASequence& c_term_ion,
			Int charge, Residue::ResidueType n_term_type, DoubleReal& n_term1, DoubleReal& c_term1, DoubleReal& n_term2, DoubleReal& c_term2*/
			// add up charges from the ions
			DoubleReal n_term_sum(0);
			for (Size i = 0; i <= n_term_ion.size(); ++i)
			{
				n_term_sum += bb_charge_[i];
				if (i != n_term_ion.size())
				{
					n_term_sum += sc_charge_[i];
				}
			}
			DoubleReal c_term_sum(0);
			for (Size i = n_term_ion.size() + 1; i != bb_charge_.size(); ++i)
			{
				c_term_sum += bb_charge_[i];
			}

			for (Size i = n_term_ion.size(); i != sc_charge_.size(); ++i)
			{
				c_term_sum += sc_charge_[i];
			}
			
			if (n_term_sum > 2)
			{
				n_term2 = 1;
				n_term1 = 0;
			}
			else
			{
				if (n_term_sum > 1)
				{
					n_term2 = n_term_sum - 1;
					n_term1 = 1 - n_term2;
				}
				else
				{
					n_term2 = 0;
					n_term1 = n_term_sum;
				}
			}

			if (c_term_sum > 2)
      {
        c_term2 = 1;
        c_term1 = 0;
      }
      else
      {
        if (c_term_sum > 1)
        {
          c_term2 = c_term_sum - 1;
          c_term1 = 1 - c_term2;
        }
				else
				{
					c_term2 = 0;
					c_term1 = c_term_sum;
				}
      }

			/*
			if (n_term_ion.size() == 2)
			{
				n_term1 /= 10.0;
				n_term2 /= 10.0;
			}
			*/

			
		}
		return;
	}

	void ProtonDistributionModel::getLeftAndRightGBValues_(const AASequence& peptide, DoubleReal& left_gb, DoubleReal& right_gb, Size position)
	{
		// TODO test if position out of range
		if (position == 0)
		{
			left_gb = (DoubleReal)param_.getValue("gb_bb_l_NH2");
			right_gb = peptide[position].getBackboneBasicityRight();
			return;
			//cerr << position << " " << left_gb << " " << right_gb << endl;
		}
		else
		{
			if (position == peptide.size())
			{
				left_gb = peptide[position - 1].getBackboneBasicityLeft();
				right_gb = (DoubleReal)param_.getValue("gb_bb_r_COOH");
				return;
				//cerr << position << " " << left_gb << " " << right_gb << endl;
			}
			else
			{
				left_gb = peptide[position - 1].getBackboneBasicityLeft();
				right_gb = peptide[position].getBackboneBasicityRight();
				return;
				//cerr << position << " " << left_gb << " " << right_gb << endl;
			}
		}
		return;
	}

} // namespace OpenMS


