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

#include <OpenMS/ANALYSIS/ID/ProtonDistributionModel.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <cmath>

#define PROTON_LOSS_FACTOR 0.3
#define PROTON_LOSS_DIFF 18000
#define PROTON_LOSS_POWER 1

using namespace std;

namespace OpenMS 
{
	ProtonDistributionModel::ProtonDistributionModel()
		:	E_(0),
			E_c_term_(0),
			E_n_term_(0)
	{
		init_();
	}

	ProtonDistributionModel::~ProtonDistributionModel()
	{
	}

	ProtonDistributionModel::ProtonDistributionModel(const ProtonDistributionModel& model)
		: sc_charge_(model.sc_charge_),
	    bb_charge_(model.bb_charge_),
			sc_charge_full_(model.sc_charge_full_),
			bb_charge_full_(model.bb_charge_full_),
			gb_sc_(model.gb_sc_),
			gb_bb_l_(model.gb_bb_l_),
			gb_bb_r_(model.gb_bb_r_),
			E_(model.E_),
			E_c_term_(model.E_c_term_),
			E_n_term_(model.E_n_term_)
	{
	}
	
	ProtonDistributionModel& ProtonDistributionModel::operator = (const ProtonDistributionModel& model)
	{
		if (this != &model)
		{
			sc_charge_ = model.sc_charge_;
      bb_charge_  = model.bb_charge_;
      sc_charge_full_ = model.sc_charge_full_;
      bb_charge_full_ = model.bb_charge_full_;
      gb_sc_ = model.gb_sc_;
      gb_bb_l_ = model.gb_bb_l_;
      gb_bb_r_ = model.gb_bb_r_;
      E_ = model.E_;
      E_c_term_ = model.E_c_term_;
      E_n_term_ = model.E_n_term_;
		}
		return *this;
	}
	
	void ProtonDistributionModel::setPeptideProtonDistribution(const HashMap<Size, double>& bb_charge, const HashMap<Size, double>& sc_charge)
	{
		bb_charge_full_ = bb_charge;
		sc_charge_full_ = sc_charge;
	}

	void ProtonDistributionModel::getProtonDistribution( HashMap<Size, double>& bb_charges,
															HashMap<Size, double>& sc_charges,
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
															
	void ProtonDistributionModel::init_()
	{
		//HashMap<String, double> gb_sc;
	  gb_sc_["D"] = 784.0;
	  gb_sc_["E"] = 790.0;
	  gb_sc_["H"] = 927.84;
	  gb_sc_["K"] = 926.74;
	  gb_sc_["M"] = 830.0;
	  gb_sc_["N"] = 864.94;
	  gb_sc_["Q"] = 865.25;
	  gb_sc_["R"] = 1000.0;
	  gb_sc_["S"] = 775.0;
	  gb_sc_["T"] = 780.0;
	  gb_sc_["W"] = 909.53;
	  gb_sc_["Y"] = 790.0;
	  //HashMap<String, double> gb_bb_l;
	  gb_bb_l_["A"] = 881.82;
	  gb_bb_l_["C"] = 880.99; // CmC
	  gb_bb_l_["D"] = 880.02;
	  gb_bb_l_["E"] = 880.10;
	  gb_bb_l_["F"] = 881.08;
	  gb_bb_l_["G"] = 881.17;
	  gb_bb_l_["H"] = 881.27;
	  gb_bb_l_["I"] = 880.99;
	  gb_bb_l_["K"] = 880.06;
	  gb_bb_l_["L"] = 881.88;
	  gb_bb_l_["M"] = 881.38;
	  gb_bb_l_["N"] = 881.18;
	  gb_bb_l_["P"] = 881.25;
	  gb_bb_l_["Q"] = 881.50;
	  gb_bb_l_["R"] = 882.98;
		gb_bb_l_["S"] = 881.08;
	  gb_bb_l_["T"] = 881.14;
	  gb_bb_l_["V"] = 881.17;
	  gb_bb_l_["W"] = 881.31;
	  gb_bb_l_["Y"] = 881.20;
	  gb_bb_l_["NH2"] = 916.84;
	  //HashMap<String, double> gb_bb_r;
	  gb_bb_r_["A"] = 0.00;
	  gb_bb_r_["C"] = 0.12; // CmC
	  gb_bb_r_["D"] = -0.63;
	  gb_bb_r_["E"] = -0.39;
	  gb_bb_r_["F"] = 0.03;
	  gb_bb_r_["G"] = 0.92;
	  gb_bb_r_["H"] = -0.19;
	  gb_bb_r_["I"] = -1.17;
	  gb_bb_r_["K"] = -0.71;
	  gb_bb_r_["L"] = -0.09;
	  gb_bb_r_["M"] = 0.30;
	  gb_bb_r_["N"] = 1.56;
	  gb_bb_r_["P"] = 11.75;
	  gb_bb_r_["Q"] = 4.10;
	  gb_bb_r_["R"] = 6.28;
	  gb_bb_r_["S"] = 0.98;
	  gb_bb_r_["T"] = 1.21;
	  gb_bb_r_["V"] = -0.90;
 	 	gb_bb_r_["W"] = 0.10;
 	 	gb_bb_r_["Y"] = -0.38;
 	 	gb_bb_r_["COOH"] = -95.82;
 	 	gb_bb_r_["b-ion"] = 36.46;
	  gb_bb_r_["a-ion"] = 46.85;
		
	}

	void ProtonDistributionModel::calculateProtonDistribution_(const AASequence& peptide, 
																								int charge, Residue::ResidueType res_type, 
																								bool fixed_proton, 
																								Size cleavage_site,
																								bool use_most_basic_site)
	{

	//cerr << "calculateProtonDistribution_(" << peptide << ", " << charge << ", " << res_type << ", " << fixed_proton << ", " << cleavage_site << ", " << use_most_basic_site << ")" << endl;

	Size most_basic_site(0);
	bool most_basic_site_sc(false);

	if (!use_most_basic_site)
	{
		bb_charge_.clear();
		sc_charge_.clear();
	}
	else
	{
		// find the most basic site
		double max_prob(0);
		for (Size i = 0; i != bb_charge_.size(); ++i)
		{
			if (bb_charge_[i] > max_prob)
			{
				max_prob = bb_charge_[i];
				most_basic_site = i;
			}
		}
		for (Size i = 0; i != sc_charge_.size(); ++i)
		{
			if (sc_charge_[i] > max_prob)
			{
				max_prob = sc_charge_[i];
				most_basic_site = i;
				most_basic_site_sc = true;
			}
		}

		bb_charge_.clear();
		sc_charge_.clear();
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

	const double T(500.0);
	
	//HashMap<Size, double> sc_charge; // side chain charges
	//HashMap<Size, double> bb_charge; // back bone charges

	for (Size i = 0; i != peptide.size(); ++i)
	{
		sc_charge_[i] = 0;
		bb_charge_[i] = 0;
	}
	bb_charge_[peptide.size()] = 0;

	// single charged
	double q(0), sum_E(0), sum_E_n_term(0), sum_E_c_term(0); // Zustandsumme
	if (charge == 1)
	{
		for (Size i = 0; i != peptide.size(); ++i)
		{
			String aa(peptide[i]->getOneLetterCode());
			
			// backbone energy
			if (i == 0)
			{
				double E = -(gb_bb_l_["NH2"] + gb_bb_r_[aa]);
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
						E = -(gb_bb_l_[aa] + gb_bb_r_["b-ion"]);
					}
					else
					{
						if (res_type == Residue::AIon)
						{
							E = -(gb_bb_l_[aa] + gb_bb_r_["a-ion"]);
						}
						else
						{
							E = -(gb_bb_l_[aa] + gb_bb_r_["COOH"]);
						}
					}
					q += exp(-E * 1000 / (Constants::R * T));
					E = -(gb_bb_l_[peptide[i - 1]->getOneLetterCode()] + gb_bb_r_[aa]);
					q += exp(-E * 1000/ (Constants::R * T));
				}
				else
				{
					// normal internal backbone position
					double E = -(gb_bb_l_[peptide[i - 1]->getOneLetterCode()] + gb_bb_r_[aa]);
					q += exp(-E * 1000 / (Constants::R * T));
				}
			}
			
			// side chains
			if (gb_sc_.has(aa))
			{
				double E = -gb_sc_[aa];
				q += exp(-E * 1000 / (Constants::R * T));
			}
		}
	
		#if 0
		cout << "Q=" << q << endl;
		#endif
		//cerr << "E: " << sum_E << endl;

		// calculate the availabilities
		for (Size i = 0; i != peptide.size(); ++i)
		{
			String aa(peptide[i]->getOneLetterCode());
			// backbone
			if (i == 0)
			{
				double E = -(gb_bb_l_["NH2"] + gb_bb_r_[aa]);
				bb_charge_[i] = exp(-E * 1000 /(Constants::R * T))/q;
				sum_E += exp(-E * 1000/Constants::R/T);
			}
			else
			{
				String aal(peptide[i-1]->getOneLetterCode());
				if (i == peptide.size() - 1)
				{
					double E(0);
					
					if (res_type == Residue::BIon)
					{
						E = -(gb_bb_l_[aa] + gb_bb_r_["b-ion"]);
					}
					else
					{
						if (res_type == Residue::AIon)
						{
							E = -(gb_bb_l_[aa] + gb_bb_r_["a-ion"]);
						}
						else
						{
							E = -(gb_bb_l_[aa] + gb_bb_r_["COOH"]);
						}
					}
					bb_charge_[i] = exp(-E * 1000/(Constants::R * T))/q;
					sum_E += exp(-E * 1000/Constants::R/T);
					E = -(gb_bb_l_[aal] + gb_bb_r_[aa]);
					bb_charge_[i+1] = exp(-E * 1000 /(Constants::R * T))/q;
					sum_E += exp(-E * 1000/Constants::R/T);
				}
				else
				{
					// normal backbone position
					double E = -(gb_bb_l_[aal] + gb_bb_r_[aa]);
					bb_charge_[i] = exp(-E * 1000 /(Constants::R * T))/q;
					sum_E += exp(-E * 1000/Constants::R/T);
				}
			}
	
			// side chains
			if (gb_sc_.has(aa))
			{
				double E = -gb_sc_[aa];
				sc_charge_[i] = exp(-E * 1000 / (Constants::R * T))/q;
				sum_E += exp(-E * 1000/Constants::R/T);
			}
		}
	}

	// fixed proton
	//
	// if two proton are available one proton is kept at the cleavage site
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
				gb_j = gb_bb_l_["NH2"] + gb_bb_r_[peptide[fixed_site]->getOneLetterCode()];
			}
			else
			{
				gb_j	= gb_bb_l_[peptide[fixed_site - 1]->getOneLetterCode()] + gb_bb_r_[peptide[fixed_site]->getOneLetterCode()];
			}
		}
		else
		{
			gb_j = gb_sc_[peptide[fixed_site]->getOneLetterCode()];
		}
		
		for (Size i = 0; i <= peptide.size(); ++i)
		{
			double gb_i(0);
      // proton 1 at N-terminus
      if (i == 0 || (i == cleavage_site && use_most_basic_site))
      {
        gb_i = gb_bb_l_["NH2"] + gb_bb_r_[peptide[i]->getOneLetterCode()];
      }
      else
      {
        String aa_i_l(peptide[i-1]->getOneLetterCode());
        // proton 1 at N-terminus
        if (i == peptide.size())
        {
					if (res_type == Residue::BIon)
					{
						gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["b-ion"];
					}
					else
					{
						if (res_type == Residue::AIon)
						{
							gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["a-ion"];
						}
						else
						{
          		gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["COOH"];
						}
					}
        }
        else
        {
          // proton 1 at backbone
          gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_[peptide[i]->getOneLetterCode()];
        }
      }

			if (!fixed_site_sc)
			{
				if (i != fixed_site)
				{
					int r_ij(abs((int)i - (int)(fixed_site)));
					q += exp(-(-gb_i - gb_j + 47.0/r_ij) * 1000 / (Constants::R * T) - 500);
	
					double gb_i_sc(0);
					if (i != peptide.size())
					{
						if (gb_sc_.has(peptide[i]->getOneLetterCode()))
						{
							gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
							q += exp(-(-gb_i_sc - gb_j + 47.0/(r_ij + 1)) * 1000 /(Constants::R * T) - 500);
						}
					}
				}
				else
				{
					// last chance: the proton i is located at side chain of cleavage site
					if (i != peptide.size())
					{
						if (gb_sc_.has(peptide[i]->getOneLetterCode()))
						{
							double gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
							q += exp(-(-gb_i - gb_i_sc + 47.0) * 1000 / (Constants::R * T) - 500);
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
				q += exp(-(-gb_i - gb_j + 47.0/(r_ij + 1)) * 1000 / (Constants::R * T) - 500);

				// only side chain site different from fixed one
				if (i != fixed_site && i != peptide.size())
				{
					double gb_i_sc(0);
					gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
					q += exp(-(-gb_i_sc - gb_j + 47.0/(r_ij + 2)) * 1000 /(Constants::R * T) - 500);
				}
			}
		}

		// calculate availablities
		for (Size i = 0; i <= peptide.size(); ++i)
		{
			double gb_i(0);
			if (i == 0 || (i == cleavage_site && use_most_basic_site))
			{
				gb_i = gb_bb_l_["NH2"] + gb_bb_r_[peptide[i]->getOneLetterCode()];
			}
			else
			{
				String aa_i_l(peptide[i - 1]->getOneLetterCode());
				if (i == peptide.size())
				{
					if (res_type == Residue::BIon)
					{
						gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["b-ion"];
					}
					else
					{
						if (res_type == Residue::AIon)
						{
							gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["a-ion"];
						}
						else
						{
							gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["COOH"];
						}
					}
				}
				else
				{
					gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_[peptide[i]->getOneLetterCode()];
				}
			}
			
			if (!fixed_site_sc)
			{
				if (i != fixed_site)
				{
					int r_ij(abs((int)i - (int)(fixed_site)));
					double prob = exp(-(-gb_i - gb_j + 47.0/r_ij) * 1000 / (Constants::R * T) - 500)/q;
					bb_charge_[i] += prob;
	
					//double add_E = exp((-gb_i - gb_j) * prob);
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
						if (gb_sc_.has(peptide[i]->getOneLetterCode()))
						{
							gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
							double prob = exp(-(-gb_i_sc - gb_j + 47.0/(r_ij + 1)) * 1000 /(Constants::R * T) - 500)/q;
							sc_charge_[i] += prob;
	
							//double add_E = exp((-gb_i_sc - gb_j) * prob);
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
						if (gb_sc_.has(peptide[i]->getOneLetterCode()))
						{
							gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
							double prob = exp(-(-gb_i_sc - gb_j + 47.0) * 1000 /(Constants::R * T) - 500)/q;
							sc_charge_[i] += prob;

							//double add_E = exp((-gb_i_sc - gb_j) * prob);
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
				double prob = exp(-(-gb_i - gb_j + 47.0/(r_ij + 1)) * 1000 / (Constants::R * T) - 500)/q;
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
					if (gb_sc_.has(peptide[i]->getOneLetterCode()))
					{
						gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
						double prob = exp(-(-gb_i_sc - gb_j + 47.0/(r_ij + 2)) * 1000 / (Constants::R * T) - 500)/q;
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
		for (Size i = 0; i <= peptide.size(); ++i)
		{
			for (Size j = i; j <= peptide.size(); ++j)
			{
				double gb_i(0), gb_j(0);
				// proton 1 at N-terminus
				if (i == 0)
				{
					gb_i = gb_bb_l_["NH2"] + gb_bb_r_[peptide[i]->getOneLetterCode()];
				}
				else
				{
					String aa_i_l(peptide[i-1]->getOneLetterCode());
					// proton 1 at N-terminus
					if (i == peptide.size())
					{
						if (res_type == Residue::BIon)
						{
							gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["b-ion"];
						}
						else
						{
							if (res_type == Residue::AIon)
							{
								gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["a-ion"];
							}
							else
							{
								gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["COOH"];
							}
						}
					}
					else
					{
						// proton 1 at backbone
						gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_[peptide[i]->getOneLetterCode()];
					}
				}
				// proton 2 at N-terminus
				if (j == 0)
				{
					gb_j = gb_bb_l_["NH2"] + gb_bb_r_[peptide[j]->getOneLetterCode()];
				}
				else
				{
					String aa_j_l(peptide[j-1]->getOneLetterCode());
					// proton 2 at C-terminus
					if (j == peptide.size())
					{
						if (res_type == Residue::BIon)
						{
							gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["b-ion"];
						}
						else
						{
							if (res_type == Residue::AIon)
							{
								gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["a-ion"];
							}
							else
							{
								gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["COOH"];
							}
						}
					}
					else
					{
						// proton 2 at backbone
						gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_[peptide[j]->getOneLetterCode()];
					}
				}
				if (i != j)
				{
					// distance of protons
					int r_ij(abs((int)i - (int)j));
					q += exp(-(-gb_i - gb_j + 47.0/r_ij) * 1000 /(Constants::R * T) - 500);
					//cerr << "1.\t" << -(-gb_i - gb_j + 47.0/r_ij) * 1000/(Constants::R * T) << endl;
					++count;

					double gb_i_sc(0), gb_j_sc(0);
					if (i != peptide.size())
					{
						// side chain of proton 1
						if (gb_sc_.has(peptide[i]->getOneLetterCode()))
						{
							gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
							q += exp(-(-gb_i_sc - gb_j + 47.0/(r_ij + 1)) * 1000 /(Constants::R * T) - 500);
							//cerr << "2.\t" << -(-gb_i_sc - gb_j + 47.0/(r_ij + 1)) * 1000/(Constants::R * T) << endl;
							++count;
						}
					}
					if (j != peptide.size())
					{
						// side chain of proton 2
						if (gb_sc_.has(peptide[j]->getOneLetterCode()))
						{
							gb_j_sc = gb_sc_[peptide[j]->getOneLetterCode()];
							q += exp(-(-gb_i - gb_j_sc + 47.0/(r_ij + 1)) * 1000 /(Constants::R * T) - 500);
							//cerr << "3.\t" << -(-gb_i - gb_j_sc + 47.0/(r_ij + 1)) * 1000 /(Constants::R * T) - 500 << endl;
							++count;
							// both at side chain?
							if (gb_i_sc != 0)
							{
								q += exp(-(-gb_i_sc - gb_j_sc + 47.0/(r_ij + 2)) * 1000/(Constants::R * T) - 500);
								//cerr << "4.\t" << -(-gb_i_sc - gb_j_sc + 47.0/(r_ij + 2)) * 1000 /(Constants::R * T) - 500 << endl;
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
						if (gb_sc_.has(peptide[i]->getOneLetterCode()))
						{
							double gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
							q += exp(-(-gb_i - gb_i_sc + 47.0) * 1000 / (Constants::R * T) - 500);
							//cerr << "5.\t" << -(-gb_i - gb_i_sc + 47.0) * 1000/ (Constants::R * T) -500 << endl;
							++count;
						}
					}
				}
			}
		}

		#if 0
		cout << "Q=" << q << ", #microstates=" << count << endl;
		#endif
		// calculate availabilities
		for (Size i = 0; i <= peptide.size(); ++i)
		{
			for (Size j = i; j <= peptide.size(); ++j)
			{
				double gb_i(0), gb_j(0);
				// calculate the backbone proton gb's
				// N-terminus
				if (i == 0)
				{
					gb_i = gb_bb_l_["NH2"] + gb_bb_r_[peptide[i]->getOneLetterCode()];
				}
				else
				{
					String aa_i_l(peptide[i-1]->getOneLetterCode());
					
					// C-terminus
					if (i == peptide.size())
					{
						if (res_type == Residue::BIon)
						{
							gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["b-ion"];
						}
						else
						{
							if (res_type == Residue::AIon)
							{
								gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["a-ion"];
							}
							else
							{
								gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_["COOH"];
							}
						}
					}
					else
					{
						// internal BB gb's
						gb_i = gb_bb_l_[aa_i_l] + gb_bb_r_[peptide[i]->getOneLetterCode()];
					}
				}
				// N-terminus
				if (j == 0)
				{
					gb_j = gb_bb_l_["NH2"] + gb_bb_r_[peptide[j]->getOneLetterCode()];
				}
				else
				{
					String aa_j_l(peptide[j-1]->getOneLetterCode());
					// C-terminus
					if (j == peptide.size())
					{
						if (res_type == Residue::BIon)
						{
							gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["b-ion"];
						}
						else
						{
							if (res_type == Residue::AIon)
							{
								gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["a-ion"];
							}
							else
							{
								gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_["COOH"];
							}
						}
					}
					else
					{
						gb_j = gb_bb_l_[aa_j_l] + gb_bb_r_[peptide[j]->getOneLetterCode()];
					}
				}

				// protons at different residues
				if (i != j)
				{
					// distance of the protons
					int r_ij(abs((int)i - (int)j));
					// calc probability
					double prob = exp(-(-gb_i - gb_j + 47.0/r_ij) * 1000 / (Constants::R * T) - 500)/q;
					// add prob to site of first proton
					bb_charge_[i] += prob;
					// add to apperent GB
					bb_charge_[j] += prob;

					// side chains
					double gb_i_sc(0), gb_j_sc(0);
					if (i != peptide.size())
					{
						if (gb_sc_.has(peptide[i]->getOneLetterCode()))
						{
							gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
							double prob = exp(-(-gb_i_sc - gb_j + 47.0/(r_ij + 1)) * 1000 /(Constants::R * T) - 500)/q;
							sc_charge_[i] += prob;
							bb_charge_[j] += prob;
						}
					}
					
					if (j != peptide.size())
					{
						if (gb_sc_.has(peptide[j]->getOneLetterCode()))
						{
							gb_j_sc = gb_sc_[peptide[j]->getOneLetterCode()];
							double prob = exp(-(-gb_i - gb_j_sc + 47.0/(r_ij + 1)) * 1000 /(Constants::R * T) - 500)/q;
							bb_charge_[i] += prob;
							sc_charge_[j] += prob;

							// both protons at sidechains
							if (gb_i_sc != 0)
							{
								double prob = exp(-(-gb_i_sc - gb_j_sc + 47.0/(r_ij + 2)) * 1000 /(Constants::R * T) - 500)/q;
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
						if (gb_sc_.has(peptide[i]->getOneLetterCode()))
						{
							double gb_i_sc = gb_sc_[peptide[i]->getOneLetterCode()];
							double prob = exp(-(-gb_i - gb_i_sc + 47.0) * 1000 / (Constants::R * T) - 500)/q;
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
			cerr << i << ".\t" << peptide[i]->getThreeLetterCode() << ": " << sc_charge_[i] << endl;
			sum += sc_charge_[i];
		}
		else
		{
			cerr << i << ".\t" << peptide[i]->getThreeLetterCode() << ": 0" << endl;
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
				n_term1 = n_term_kapp / (n_term_kapp + c_term_kapp);
				c_term1 = c_term_kapp / (n_term_kapp + c_term_kapp);

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
	
				// calc proton distribution of second proton with other one at most basic site fixed
				calculateProtonDistribution_(peptide, 2, Residue::Full, false, n_term_ion.size(), true);

				double singly_charged(0);
				for (Size i = 0; i != n_term_ion.size(); ++i)
				{
					n_term2 += bb_charge_[i] * p_n;
					singly_charged += bb_charge_[i] * p_c;
					if (sc_charge_.has(i))
					{
						n_term2 += sc_charge_[i] * p_n;
						singly_charged += sc_charge_[i] * p_c;
					}
				}

				for (Size i = n_term_ion.size(); i <= peptide.size(); ++i)
				{
					c_term2 += bb_charge_[i] * p_c;
					singly_charged += bb_charge_[i] * p_n;
					if (sc_charge_.has(i - 1))
					{
						c_term2 += sc_charge_[i - 1] * p_c;
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

				//cerr << "charge=2, n_term1=" << n_term1 << ", n_term2=" << n_term2 << ", c_term1=" << c_term1 << ", c_term2=" << c_term2 << endl;
			}
			else
			{
				if (type == ChargeRemote || type == SideChain)
				{
					// TODO ranges correct? Missing some sites of the peptide
					double n_term_sum(0), c_term_sum(0);
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

					double sum(0);
					sum += n_term1 + n_term2 + c_term1 + c_term2;
					n_term1 /= sum;
					n_term2 /= sum;
					c_term1 /= sum;
					c_term2 /= sum;
						
				}
				else
				{
					cerr << "calcChargeStateIntensities_: unknown fragmentation type (" << type << ")" << endl;
				}
			}
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
			for (Size i = peptide.size() - ion.size(); i != peptide.size(); ++i)
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
			for (Size i = 0; i <= ion.size(); ++i)
			{
				sum += bb_charge_[i];
				if (sc_charge_.has(i))
				{
					sum += sc_charge_[i];
				}
			}
		}
		
		//cout << "charge: " << sum << endl;
		
		
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
			//cout << "charge ratio: " << tmp/(tmp+1) << " " << 1/(tmp+1) << " " << tmp << " " << pa << endl;
		}
		
		return ints;
	}

} // namespace OpenMS

