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

#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <cmath>
#include <numeric>
#include <cstdlib>

#define PROTON_LOSS_FACTOR 0.3
#define PROTON_LOSS_DIFF 18000
#define PROTON_LOSS_POWER 1
#define COULOMB_REPULSION (double)47.0  // from zhang: 47.0 kJ/mol


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
		//defaults_.setValue("temperature", 500.0, "Temperature term ", true);

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
	
	void ProtonDistributionModel::setPeptideProtonDistribution(const Map<UInt, double>& bb_charge, const Map<UInt, double>& sc_charge)
	{
		bb_charge_full_ = bb_charge;
		sc_charge_full_ = sc_charge;
	}

	void ProtonDistributionModel::getProtonDistribution( Map<UInt, double>& bb_charges,
															Map<UInt, double>& sc_charges,
															const AASequence& peptide,
															int charge,
															Residue::ResidueType res_type)
	{
		calculateProtonDistribution_(peptide, charge, res_type);
		bb_charges = bb_charge_;
		sc_charges = sc_charge_;
	}

	void ProtonDistributionModel::getChargeStateIntensities(const AASequence& peptide,
													                                const AASequence& n_term_ion,
																													const AASequence& c_term_ion,
																													int charge,
																													Residue::ResidueType n_term_type,
																													double& n_term1,
																													double& c_term1,
																													double& n_term2,
																													double& c_term2,
																													FragmentationType type)
	{
		calcChargeStateIntensities_(peptide, n_term_ion, c_term_ion, charge, n_term_type, n_term1, c_term1, n_term2, c_term2, type);
	}
	
	void ProtonDistributionModel::calculateProtonDistribution_(const AASequence& peptide, 
																								int charge, Residue::ResidueType res_type, 
																								bool fixed_proton, 
																								UInt cleavage_site,
																								bool use_most_basic_site)
	{
		//cerr << "calculateProtonDistribution_(" << peptide << ", charge=" << charge	<< ", fixed_proton=" << fixed_proton << ", cleavage_site=" << cleavage_site << ", use_most_basic_site=" << use_most_basic_site << ")" << endl;
		if (charge > 2)
		{
			vector<double> gb_bb, gb_sc;
			double gb_left_n_term(0), gb_right_n_term(0);
			getLeftAndRightGBValues_(peptide, gb_left_n_term, gb_right_n_term, 0);
			gb_bb.push_back(gb_left_n_term + gb_right_n_term);
			UInt count(1);
			for (AASequence::ConstIterator it = peptide.begin(); it != peptide.end(); ++it, ++count)
			{
				double gb(0), gb_left(0), gb_right(0);
				getLeftAndRightGBValues_(peptide, gb_left, gb_right, count);
				
				gb = gb_left + gb_right;
				gb_bb.push_back(gb);

				gb_sc.push_back(it->getSideChainBasicity());
			}
			
			// now distribute the charges until no site has more than 1.0 proton
			vector<double> bb_coulomb(peptide.size() + 1, 0.0), sc_coulomb(peptide.size(), 0.0);
			Int actual_charge(charge);
			set<UInt> sc_sites, bb_sites;
			while (true)
			{
				//cerr << "#proton remaining: " << actual_charge << endl;
				vector<double> k_bb(peptide.size() + 1, 0.0), k_sc(peptide.size(), 0.0);
				count = 0;
				double sum_k(0);
				for (vector<double>::const_iterator it = gb_bb.begin(); it != gb_bb.end(); ++it, ++count)
				{
					if (bb_sites.find(count) == bb_sites.end())
					{
						k_bb[count] = exp((*it - bb_coulomb[count]) * 1000.0 / Constants::R / 500.0);
						sum_k += k_bb[count];
						//cerr << k_bb[count] << endl;
					}
				}

				count = 0;
				for (vector<double>::const_iterator it = gb_sc.begin(); it != gb_sc.end(); ++it, ++count)
				{
					if (sc_sites.find(count) == sc_sites.end())
					{
						k_sc[count] = exp((*it - sc_coulomb[count]) * 1000.0 / Constants::R / 500.0);
						sum_k += k_sc[count];
						//cerr << k_sc[count] << endl;
					}
				}

				//cerr << "sum_k: " << sum_k << endl;

				vector<double> p_bb(peptide.size() + 1, 1.0), p_sc(peptide.size(), 1.0);
				count = 0;
				for (vector<double>::const_iterator it = k_bb.begin(); it != k_bb.end(); ++it, ++count)
				{
					if (bb_sites.find(count) == bb_sites.end())
					{
						p_bb[count] = (double)actual_charge * *it / sum_k;
						bb_charge_[count] = p_bb[count];
					}
					//cerr << "BB" << count << ": " << p_bb[count] << endl;
				}

				count = 0;
				for (vector<double>::const_iterator it = k_sc.begin(); it != k_sc.end(); ++it, ++count)
				{
					if (sc_sites.find(count) == sc_sites.end())
					{
						p_sc[count] = (double)actual_charge * *it / sum_k;
						sc_charge_[count] = p_sc[count];
					}
					//cerr << "SC" << count << ": " << p_sc[count] << endl;
				}
		
				// check if there is a site containing more than one proton
				for (UInt i = 0; i != p_bb.size(); ++i)
				{
					if (p_bb[i] > 1.0)
					{
						bb_charge_[i] = 1.0;
						bb_sites.insert(i);
						//cerr << "BackbonePosition " << i << " has charge > 1" << endl;
						--actual_charge;
					}
				}
				for (UInt i = 0; i != p_sc.size(); ++i)
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
				for (UInt i = 0; i != gb_bb.size(); ++i)
				{
					// check if the site is not occupied by a "complete" proton
					if (bb_sites.find(i) == bb_sites.end())
					{
						double coulomb_sum(0);
						for (set<UInt>::const_iterator it = bb_sites.begin(); it != bb_sites.end(); ++it)
						{
							// calculate the distance between occupied site and this backbone site
							UInt pos = *it;
							UInt diff = (pos > i) ? pos - i : i - pos;
							coulomb_sum += COULOMB_REPULSION / (double)diff;
						}

						for (set<UInt>::const_iterator it = sc_sites.begin(); it != sc_sites.end(); ++it)
						{
							// calculate the distance between occupied side chain and this backbone site
							UInt pos = *it;
							UInt diff = (pos > i) ? pos -i : i - pos;
							++diff; // bond to the side chain counts extra
							coulomb_sum += COULOMB_REPULSION / (double)diff;
						}
						bb_coulomb[i] = coulomb_sum;
						//cerr << "BB coulomb" << i << ": " << coulomb_sum << endl;
					}
				}
				for (UInt i = 0; i != gb_sc.size(); ++i)
				{
					if (sc_sites.find(i) == sc_sites.end())
					{
						double coulomb_sum(0);
						for (set<UInt>::const_iterator it = bb_sites.begin(); it != bb_sites.end(); ++it)
						{
							UInt pos = *it;
							UInt diff = (pos > i) ? pos - i : i - pos;
							++diff;
							coulomb_sum += COULOMB_REPULSION / (double)diff;
						}
						for (set<UInt>::const_iterator it = sc_sites.begin(); it != sc_sites.end(); ++it)
						{
							UInt pos = *it;
							UInt diff = (pos > i) ? pos - i : i - pos;
							diff += 2;
							coulomb_sum += COULOMB_REPULSION / (double)diff;
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
				for (vector<double>::const_iterator it = p_bb.begin(); it != p_bb.end(); ++it)
				{
					if (*it > 1.0)
					{
						has_greater_one = true;
					}
				}

				for (vector<double>::const_iterator it = p_sc.begin(); it != p_sc.end(); ++it)
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
		
		UInt most_basic_site(0);
		bool most_basic_site_sc(false);

    double gb_bb_l_NH2 = param_.getValue("gb_bb_l_NH2");
		double gb_bb_r_COOH = param_.getValue("gb_bb_r_COOH");
		double gb_bb_r_bion = param_.getValue("gb_bb_r_b-ion");
		double gb_bb_r_aion = param_.getValue("gb_bb_r_a-ion");

		if (!use_most_basic_site)
		{
			bb_charge_.clear();
			sc_charge_.clear();
		}
		else
		{
			// find the most basic site
			double max_prob(0);
			//cerr << "bb: ";
			for (UInt i = 0; i != bb_charge_.size(); ++i)
			{
				//cerr << i << ". " << bb_charge_[i] << "; " << endl;
				if (bb_charge_[i] > max_prob)
				{
					max_prob = bb_charge_[i];
					most_basic_site = i;
				}
			}

			//cerr << endl << "sc: ";
			for (UInt i = 0; i != sc_charge_.size(); ++i)
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
			

			bb_charge_.clear();
			sc_charge_.clear();
		}

		UInt fixed_site(0);
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

		const double T(500.0);
	
		for (UInt i = 0; i != sc_charge_.size(); ++i)
		{
			sc_charge_[i] = 0;
		}

		for (UInt i = 0; i != bb_charge_.size(); ++i)
		{
			bb_charge_[i] = 0;
		}
		//bb_charge_[peptide.size()] = 0;

		// single charged
		double q(0), sum_E(0), sum_E_n_term(0), sum_E_c_term(0); // Zustandsumme
		if (charge == 1)
	{
		for (UInt i = 0; i != peptide.size(); ++i)
		{
			//String aa(peptide[i]->getOneLetterCode());
			
			// backbone energy
			if (i == 0)
			{
				//double E = -(gb_bb_l_["NH2"] + gb_bb_r_[aa]);
				double E = -(gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight());
				q += exp(-E * 1000 / (Constants::R * T));
			}
			else
			{
				if (i == peptide.size() - 1)
				{
					// position at the C-terminal end of the ion
					double E(0);
					if (res_type == Residue::BIon)
					{
						//E = -(gb_bb_l_[aa] + gb_bb_r_["b-ion"]);
						E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_bion);
					}
					else
					{
						if (res_type == Residue::AIon)
						{
							//E = -(gb_bb_l_[aa] + gb_bb_r_["a-ion"]);
							E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_aion);
						}
						else
						{
							//E = -(gb_bb_l_[aa] + gb_bb_r_["COOH"]);
							E = -(peptide[i].getBackboneBasicityLeft() + gb_bb_r_COOH);
						}
					}
					q += exp(-E * 1000 / (Constants::R * T));
					//E = -(gb_bb_l_[peptide[i - 1].getOneLetterCode()] + gb_bb_r_[aa]);
					E = -(peptide[i - 1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
					q += exp(-E * 1000/ (Constants::R * T));
				}
				else
				{
					// normal internal backbone position
					//double E = -(gb_bb_l_[peptide[i - 1].getOneLetterCode()] + gb_bb_r_[aa]);
					double E = -(peptide[i - 1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
					q += exp(-E * 1000 / (Constants::R * T));
				}
			}
			
			// side chains
			if (/*gb_sc_.has(aa)*/peptide[i].getSideChainBasicity() != 0)
			{
				//double E = -gb_sc_[aa];
				double E = -peptide[i].getSideChainBasicity();
				q += exp(-E * 1000 / (Constants::R * T));
			}
		}
	
		#if 0
		cout << "Q=" << q << endl;
		#endif
		//cerr << "E: " << sum_E << endl;

		// calculate the availabilities
		for (UInt i = 0; i != peptide.size(); ++i)
		{
			// backbone
			if (i == 0)
			{
				//double E = -(gb_bb_l_["NH2"] + gb_bb_r_[aa]);
				double E = -(gb_bb_l_NH2 + peptide[i].getBackboneBasicityRight());
				bb_charge_[i] = exp(-E * 1000 /(Constants::R * T))/q;
				sum_E += exp(-E * 1000/Constants::R/T);
			}
			else
			{
				if (i == peptide.size() - 1)
				{
					double E(0);
					
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
					bb_charge_[i] = exp(-E * 1000/(Constants::R * T))/q;
					sum_E += exp(-E * 1000/Constants::R/T);
					E = -(peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
					bb_charge_[i+1] = exp(-E * 1000 /(Constants::R * T))/q;
					sum_E += exp(-E * 1000/Constants::R/T);
				}
				else
				{
					// normal backbone position
					double E = -(peptide[i-1].getBackboneBasicityLeft() + peptide[i].getBackboneBasicityRight());
					bb_charge_[i] = exp(-E * 1000 /(Constants::R * T))/q;
					sum_E += exp(-E * 1000/Constants::R/T);
				}
			}
	
			// side chains
			if (peptide[i].getSideChainBasicity() != 0)
			{
				double E = -peptide[i].getSideChainBasicity();
				sc_charge_[i] = exp(-E * 1000 / (Constants::R * T))/q;
				sum_E += exp(-E * 1000/Constants::R/T);
			}
		}
	}

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
		
		double gb_j(0);
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
		
		for (UInt i = 0; i <= peptide.size(); ++i)
		{
			double gb_i(0);
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
					int r_ij(abs((int)i - (int)(fixed_site)));
					q += exp(-(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000 / (Constants::R * T) - 500);
	
					double gb_i_sc(0);
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
							double gb_i_sc = peptide[i].getSideChainBasicity();
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
				int r_ij = abs((int)i - (int)fixed_site);
				q += exp(-(-gb_i - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 / (Constants::R * T) - 500);

				// only side chain site different from fixed one
				if (i != fixed_site && i != peptide.size())
				{
					double gb_i_sc(0);
					gb_i_sc = peptide[i].getSideChainBasicity();
					q += exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 2)) * 1000 /(Constants::R * T) - 500);
				}
			}
		}

		// calculate availablities
		for (UInt i = 0; i <= peptide.size(); ++i)
		{
			double gb_i(0);
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
					int r_ij(abs((int)i - (int)(fixed_site)));
					double prob = exp(-(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000 / (Constants::R * T) - 500)/q;
					bb_charge_[i] += prob;
	
					double add_E = exp(gb_i * 1000 / Constants::R / T);
					if (i < fixed_site - 1)
					{
						sum_E_n_term += add_E;
					}
					else
					{
						sum_E_c_term += add_E;
					}

					double gb_i_sc(0);
					if (i != peptide.size())
					{
						if (peptide[i].getSideChainBasicity() != 0)
						{
							gb_i_sc = peptide[i].getSideChainBasicity();
							double prob = exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500)/q;
							sc_charge_[i] += prob;
	
							double add_E = exp(gb_i_sc * 1000 / Constants::R / T);
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
						double gb_i_sc(0);
						if (peptide[i].getSideChainBasicity() != 0)
						{
							gb_i_sc = peptide[i].getSideChainBasicity();
							double prob = exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION) * 1000 /(Constants::R * T) - 500)/q;
							sc_charge_[i] += prob;

							double add_E = exp(gb_i_sc * 1000 / Constants::R / T);
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
				int r_ij = abs((int)i - (int)fixed_site);
				double prob = exp(-(-gb_i - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 / (Constants::R * T) - 500)/q;
				bb_charge_[i] += prob;

				double add_E = exp(gb_i * 1000 / Constants::R / T);
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
					double gb_i_sc(0);
					if (peptide[i].getSideChainBasicity() != 0)
					{
						gb_i_sc = peptide[i].getSideChainBasicity();
						double prob = exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 2)) * 1000 / (Constants::R * T) - 500)/q;
						sc_charge_[i] += prob;

						double add_E = exp(gb_i_sc * 1000 / Constants::R / T);
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

	
	// double charged
	if (charge == 2 && !fixed_proton && !use_most_basic_site)
	{
		// calculate sum
		int count(0);
		for (UInt i = 0; i <= peptide.size(); ++i)
		{
			for (UInt j = i; j <= peptide.size(); ++j)
			{
				double gb_i(0), gb_j(0);
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
					double gb_j_l = peptide[j-1].getBackboneBasicityLeft();
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
					int r_ij(abs((int)i - (int)j));
					q += exp(-(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000 /(Constants::R * T) - 500);
					//cerr << "1.\t" << -(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000/(Constants::R * T) << endl;
					++count;

					double gb_i_sc(0), gb_j_sc(0);
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
							//double gb_i_sc = gb_sc_[peptide[i].getOneLetterCode()];
							double gb_i_sc = peptide[i].getSideChainBasicity();
							q += exp(-(-gb_i - gb_i_sc + COULOMB_REPULSION) * 1000 / (Constants::R * T) - 500);
							//cerr << "5.\t" << -(-gb_i - gb_i_sc + COULOMB_REPULSION) * 1000/ (Constants::R * T) -500 << endl;
							++count;
						}
					}
				}
			}
		}

		// calculate availabilities
		for (UInt i = 0; i <= peptide.size(); ++i)
		{
			for (UInt j = i; j <= peptide.size(); ++j)
			{
				double gb_i(0), gb_j(0);
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
					int r_ij(abs((int)i - (int)j));
					// calc probability
					double prob = exp(-(-gb_i - gb_j + COULOMB_REPULSION/r_ij) * 1000 / (Constants::R * T) - 500)/q;
					// add prob to site of first proton
					bb_charge_[i] += prob;
					// add to apperent GB
					bb_charge_[j] += prob;

					// side chains
					double gb_i_sc(0), gb_j_sc(0);
					if (i != peptide.size())
					{
						if (/*gb_sc_.has(peptide[i].getOneLetterCode())*/peptide[i].getSideChainBasicity() != 0)
						{
							//gb_i_sc = gb_sc_[peptide[i].getOneLetterCode()];
							gb_i_sc = peptide[i].getSideChainBasicity();
							double prob = exp(-(-gb_i_sc - gb_j + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500)/q;
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
							double prob = exp(-(-gb_i - gb_j_sc + COULOMB_REPULSION/(r_ij + 1)) * 1000 /(Constants::R * T) - 500)/q;
							bb_charge_[i] += prob;
							sc_charge_[j] += prob;

							// both protons at sidechains
							if (gb_i_sc != 0)
							{
								double prob = exp(-(-gb_i_sc - gb_j_sc + COULOMB_REPULSION/(r_ij + 2)) * 1000 /(Constants::R * T) - 500)/q;
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
							//double gb_i_sc = gb_sc_[peptide[i].getOneLetterCode()];
							double gb_i_sc = peptide[i].getSideChainBasicity();
							double prob = exp(-(-gb_i - gb_i_sc + COULOMB_REPULSION) * 1000 / (Constants::R * T) - 500)/q;
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
	double sum(0);
	for (unsigned int i = 0; i != peptide.size(); ++i)
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
	for (unsigned int i = 0; i  <= peptide.size(); ++i)
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

	double ProtonDistributionModel::getProtonAffinity_(const AASequence& peptide, int charge, Residue::ResidueType res_type)
	{
		//const double T(500.0);

		double pa(0);
		calculateProtonDistribution_(peptide, charge, res_type);

		//pa = Constants::R * T * log(E_);

		// new test
		double sum(0);
		if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
		{
			for (UInt i = 1; i <= peptide.size(); ++i)
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
				for (UInt i = bb_charge_.size() - peptide.size() - 1; i != bb_charge_.size(); ++i)
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

	void ProtonDistributionModel::calcChargeStateIntensities_( const AASequence& peptide, 
																									const AASequence& n_term_ion, 
																									const AASequence& c_term_ion,
																									int charge, 
																									Residue::ResidueType n_term_type,	
																									double& n_term1, 
																									double& c_term1, 
																									double& n_term2, 
																									double& c_term2,
																									FragmentationType type)
	{
	
		double n_term_kapp(0), c_term_kapp(0);
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
				//n_term1 = n_term_kapp / (n_term_kapp + c_term_kapp);
				//c_term1 = c_term_kapp / (n_term_kapp + c_term_kapp);
				//
				
				double pa_n = log(n_term_kapp);
				double pa_c = log(c_term_kapp);

				double ratio_bx_yz = exp(pa_n - pa_c);

				n_term1 = ratio_bx_yz / (1.0 + ratio_bx_yz);
				c_term1 = 1.0 / (ratio_bx_yz + 1.0);

				// of course ++ ions are not available
				n_term2 = 0;
				c_term2 = 0;
				

				//cerr << n_term_kapp << " " << c_term_kapp << " " << n_term1 << " " << c_term1 << endl;
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
				double p_n(0), p_c(0);

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
				for (UInt i = 0; i != bb_charge_.size(); ++i)
				{
					cerr << "; " << i << ". " << bb_charge_[i];
				}
				cerr << "\nSC: ";
				for (UInt i = 0; i != sc_charge_.size(); ++i)
				{
					cerr << "; " << i << ". " << sc_charge_[i];
				}
				cerr << endl;
#endif
				
				double singly_charged(0);
				for (UInt i = 0; i != n_term_ion.size(); ++i)
				{
					n_term2 += bb_charge_[i] * p_n;
					singly_charged += bb_charge_[i] * p_c;
					if (sc_charge_.has(i))
					{
						n_term2 += sc_charge_[i] * p_n;
						singly_charged += sc_charge_[i]  * p_c;
					}
				}

				for (UInt i = n_term_ion.size(); i <= peptide.size(); ++i)
				{
					c_term2 += bb_charge_[i] * p_c;
					singly_charged += bb_charge_[i] * p_n;
					if (sc_charge_.has(i - 1))
					{
						c_term2 += sc_charge_[i - 1]  * p_c;
						singly_charged += sc_charge_[i - 1] * p_n;
					}
				}

				n_term1 = singly_charged;
				c_term1 = singly_charged;

				//cerr << E_n_term_ << "\t" << E_c_term_ << "\t" << p_n << "\t" << p_c << "\t" << endl;

				// TODO normalization correct?
				double sum(0);
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
					double n_term_sum(0), c_term_sum(0);
					for (UInt i = 0; i != n_term_ion.size(); ++i)
					{
						n_term_sum += bb_charge_full_[i];
						n_term_sum += sc_charge_full_[i];
					}
					for (UInt i = n_term_ion.size(); i != peptide.size(); ++i)
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

					double sum(0);
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
			int charge, Residue::ResidueType n_term_type, double& n_term1, double& c_term1, double& n_term2, double& c_term2*/
			// add up charges from the ions
			double n_term_sum(0);
			for (UInt i = 0; i <= n_term_ion.size(); ++i)
			{
				n_term_sum += bb_charge_[i];
				if (i != n_term_ion.size())
				{
					n_term_sum += sc_charge_[i];
				}
			}
			double c_term_sum(0);
			for (UInt i = n_term_ion.size() + 1; i != bb_charge_.size(); ++i)
			{
				c_term_sum += bb_charge_[i];
			}

			for (UInt i = n_term_ion.size(); i != sc_charge_.size(); ++i)
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

	vector<double> ProtonDistributionModel::getChargeStateIntensities_(const AASequence& peptide, const AASequence& ion, int charge, Residue::ResidueType res_type)
	{
		const double T(8.314472);
		
		vector<double> ints;
		calculateProtonDistribution_(peptide, charge, res_type);
		double sum(0);
		if (res_type == Residue::YIon || res_type == Residue::XIon || res_type == Residue::ZIon)
		{
			for (UInt i = peptide.size() - ion.size(); i != peptide.size(); ++i)
			{
				sum += bb_charge_[i+1];
				if (sc_charge_.has(i))
				{
					sum += sc_charge_[i];
				}
			}
		}
		else
		{
			for (UInt i = 0; i <= ion.size(); ++i)
			{
				sum += bb_charge_[i];
				if (sc_charge_.has(i))
				{
					sum += sc_charge_[i];
				}
			}
		}
		
		if (sum < 1)
		{
			ints.push_back(1.0);
		}
		else
		{
			//double pa = getProtonAffinity_(ion, 1, res_type);
			calculateProtonDistribution_(ion, 1, res_type);
			double pa = Constants::R * T * log(E_);
			
		
			double tmp = PROTON_LOSS_FACTOR * exp(-((pow(pa - PROTON_LOSS_DIFF, PROTON_LOSS_POWER)))/Constants::R/T);
			ints.push_back(tmp/(tmp+1));
			ints.push_back(1/(tmp+1));
		}
		
		return ints;
	}

	void ProtonDistributionModel::getLeftAndRightGBValues_(const AASequence& peptide, double& left_gb, double& right_gb, UInt position)
	{
		// TODO test if position out of range
		if (position == 0)
		{
			left_gb = (double)param_.getValue("gb_bb_l_NH2");
			right_gb = peptide[position].getBackboneBasicityRight();
			//cerr << position << " " << left_gb << " " << right_gb << endl;
		}
		else
		{
			if (position == peptide.size())
			{
				left_gb = peptide[position - 1].getBackboneBasicityLeft();
				right_gb = (double)param_.getValue("gb_bb_r_COOH");
				//cerr << position << " " << left_gb << " " << right_gb << endl;
			}
			else
			{
				left_gb = peptide[position - 1].getBackboneBasicityLeft();
				right_gb = peptide[position].getBackboneBasicityRight();
				//cerr << position << " " << left_gb << " " << right_gb << endl;
			}
		}
	}

} // namespace OpenMS


