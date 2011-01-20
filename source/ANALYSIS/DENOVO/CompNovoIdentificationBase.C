// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationBase.h>

#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoring.h>

#define MIN_DOUBLE_MZ 900.0

//#define SELECT_PIVOT_DEBUG

using namespace std;

namespace OpenMS
{
	CompNovoIdentificationBase::CompNovoIdentificationBase()
		:	DefaultParamHandler("CompNovoIdentificationBase"),
		  max_number_aa_per_decomp_(0),
      tryptic_only_(true),
     	fragment_mass_tolerance_(0),
      max_number_pivot_(0),
      decomp_weights_precision_(0),
      max_mz_(2000.0),
      min_mz_(200.0),
      max_decomp_weight_(450.0),
      max_subscore_number_(30),
			max_isotope_(3)
	{
		defaults_.setValue("max_number_aa_per_decomp", 4, "maximal amino acid frequency per decomposition", StringList::create("advanced"));
		defaults_.setValue("tryptic_only", "true", "if set to true only tryptic peptides are reported");
		defaults_.setValue("precursor_mass_tolerance", 1.5, "precursor mass tolerance");
		defaults_.setValue("fragment_mass_tolerance", 0.3, "fragment mass tolerance");
		defaults_.setValue("max_number_pivot", 9, "maximal number of pivot ions to be used", StringList::create("advanced"));
		defaults_.setValue("max_subscore_number", 40, "maximal number of solutions of a subsegment that are kept", StringList::create("advanced"));
		defaults_.setValue("decomp_weights_precision", 0.01, "precision used to calculate the decompositions, this only affects cache usage!", StringList::create("advanced"));
		defaults_.setValue("double_charged_iso_threshold", 0.6, "minimal isotope intensity correlation of doubly charged ions to be used to score the single scored ions", StringList::create("advanced"));
		defaults_.setValue("max_mz", 2000.0, "maximal m/z value used to calculate isotope distributions");
		defaults_.setValue("min_mz", 200.0, "minimal m/z value used to calculate the isotope distributions");
		defaults_.setValue("max_isotope_to_score", 3, "max isotope peak to be considered in the scoring", StringList::create("advanced"));
		defaults_.setValue("max_decomp_weight", 450.0, "maximal m/z difference used to calculate the decompositions", StringList::create("advanced"));
		defaults_.setValue("max_isotope", 3, "max isotope used in the theoretical spectra to score", StringList::create("advanced"));
		defaults_.setValue("missed_cleavages", 1, "maximal number of missed cleavages allowed per peptide");
		defaults_.setValue("number_of_hits", 100, "maximal number of hits which are reported per spectrum");
		defaults_.setValue("estimate_precursor_mz", "true", "If set to true, the precursor charge will be estimated, e.g. from the precursor peaks of the ETD spectrum.\n"
																												 "The input is believed otherwise.");
		defaults_.setValidStrings("estimate_precursor_mz", StringList::create("true,false"));
		defaults_.setValue("number_of_prescoring_hits", 250, "how many sequences are kept after first rough scoring for better scoring", StringList::create("advanced"));

		// set all known modifications as restriction
		vector<String> all_mods;
		ModificationsDB::getInstance()->getAllSearchModifications(all_mods);

		defaults_.setValue("fixed_modifications", StringList::create(""), "fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'");
		defaults_.setValidStrings("fixed_modifications", all_mods);
		
		defaults_.setValue("variable_modifications", StringList::create(""), "variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'");
		defaults_.setValidStrings("variable_modifications", all_mods);
		
		defaults_.setValue("residue_set", "Natural19WithoutI", "The predefined amino acid set that should be used, see doc of ResidueDB for possible residue sets", StringList::create("advanced"));

		defaultsToParam_();
	}

	CompNovoIdentificationBase::CompNovoIdentificationBase(const CompNovoIdentificationBase& rhs)
		: DefaultParamHandler(rhs)
	{
		updateMembers_();
	}

	CompNovoIdentificationBase& CompNovoIdentificationBase::operator = (const CompNovoIdentificationBase& rhs)
	{
		if (this != &rhs)
		{
			DefaultParamHandler::operator=(rhs);
			updateMembers_();
			// TODO
		}
		return *this;
	}

	CompNovoIdentificationBase::~CompNovoIdentificationBase()
	{
	}

	void CompNovoIdentificationBase::getCIDSpectrumLight_(PeakSpectrum& spec, const String& sequence, DoubleReal prefix, DoubleReal suffix)
	{
		static DoubleReal h2o_mass = EmpiricalFormula("H2O").getMonoWeight();
  	Peak1D p;
  	DoubleReal b_pos(0.0 + prefix);
	  DoubleReal y_pos(h2o_mass + suffix);
  	for (Size i = 0; i != sequence.size() - 1; ++i)
	  {
  	  char aa(sequence[i]);
    	b_pos += aa_to_weight_[aa];

    	char aa2(sequence[sequence.size() - i - 1]);
    	y_pos += aa_to_weight_[aa2];

    	if (b_pos > min_mz_ && b_pos < max_mz_)
    	{
      	p.setPosition(b_pos + Constants::PROTON_MASS_U);
      	p.setIntensity(1.0f);
      	spec.push_back(p);
    	}

    	if (y_pos > min_mz_ && y_pos < max_mz_)
    	{
      	p.setPosition(y_pos + Constants::PROTON_MASS_U);
      	p.setIntensity(1.0f);
      	spec.push_back(p);
    	}
  	}

  	spec.sortByPosition();
		return;
	}

	void CompNovoIdentificationBase::getCIDSpectrum_(PeakSpectrum& spec, const String& sequence, Size charge, DoubleReal prefix, DoubleReal suffix)
	{
		static DoubleReal h2o_mass = EmpiricalFormula("H2O").getMonoWeight();
		static DoubleReal nh3_mass = EmpiricalFormula("NH3").getMonoWeight();
		static DoubleReal co_mass = EmpiricalFormula("CO").getMonoWeight();
		Peak1D p;
    DoubleReal b_pos(0 + prefix);
    DoubleReal y_pos(h2o_mass + suffix);
    bool b_H2O_loss(false), b_NH3_loss(false), y_NH3_loss(false);

    for (Size i = 0; i != sequence.size() - 1; ++i)
    {
      char aa(sequence[i]);
      b_pos += aa_to_weight_[aa];

      char aa2(sequence[sequence.size() - i - 1]);
      y_pos += aa_to_weight_[aa2];
      for (Size z = 1; z <= charge && z < 3; ++z)
      {
        // b-ions
				if (b_pos >= min_mz_ && b_pos <= max_mz_)
				{
					for (Size j = 0; j != max_isotope_; ++j)
					{
						if (z == 1 /*|| b_pos > MIN_DOUBLE_MZ*/)
						{
        			p.setPosition((b_pos + (DoubleReal)z * Constants::PROTON_MASS_U + (DoubleReal)j + Constants::NEUTRON_MASS_U)/(DoubleReal)z);
        			p.setIntensity(isotope_distributions_[(Size)b_pos][j] * 0.8 / (z*z));
							spec.push_back(p);
						}
					}
				}
			
        // b-ion losses
				if (b_pos - h2o_mass > min_mz_ && b_pos - h2o_mass < max_mz_)
				{	
        	if (b_H2O_loss || aa == 'S' || aa == 'T' || aa == 'E' || aa == 'D')
        	{
          	b_H2O_loss = true;
          	p.setPosition((b_pos + z * Constants::PROTON_MASS_U - h2o_mass)/z);
          	p.setIntensity(0.02 / (DoubleReal)(z*z));
						if (z == 1/* || b_pos > MIN_DOUBLE_MZ*/)
						{
          		spec.push_back(p);
						}
        	}
        	if (b_NH3_loss || aa == 'Q' || aa == 'N' || aa == 'R' || aa == 'K')
        	{
          	b_NH3_loss = true;
          	p.setPosition((b_pos + z * Constants::PROTON_MASS_U - nh3_mass)/z);
          	p.setIntensity(0.02 / (DoubleReal)(z*z));

						if (z == 1/* || b_pos > MIN_DOUBLE_MZ*/)
						{
         			spec.push_back(p);
						}
        	}
				}

				// a-ions only for charge 1
				if (z == 1)
				{
					if (b_pos - co_mass > min_mz_ && b_pos - co_mass < max_mz_)
					{
        		// a-ions
        		p.setPosition((b_pos + z * Constants::PROTON_MASS_U - co_mass)/(DoubleReal)z);
       			p.setIntensity(0.1f);
        		spec.push_back(p);
					}
				}
			

				
				if (y_pos > min_mz_ && y_pos < max_mz_)
				{
        	// y-ions
					for (Size j = 0; j != max_isotope_; ++j)
					{
						if (z == 1/* || y_pos > MIN_DOUBLE_MZ*/)
						{
							p.setPosition((y_pos + (DoubleReal)z * Constants::PROTON_MASS_U + (DoubleReal)j * Constants::NEUTRON_MASS_U)/(DoubleReal)z);
        			p.setIntensity(isotope_distributions_[(Size)y_pos][j] /(DoubleReal) (z*z));
        			spec.push_back(p);
						}
					}
		
					// H2O loss
          p.setPosition((y_pos + z * Constants::PROTON_MASS_U - h2o_mass)/(DoubleReal)z);
          p.setIntensity(0.1 / (DoubleReal)(z*z));
					if (aa2 == 'Q') // pyroglutamic acid formation
					{
						p.setIntensity(0.5f);
					}	
					if (z == 1/* || y_pos > MIN_DOUBLE_MZ*/)
					{
					   spec.push_back(p);
					}

					// NH3 loss
        	if (y_NH3_loss || aa2 == 'Q' || aa2 == 'N' || aa2 == 'R' || aa2 == 'K')
        	{
          	y_NH3_loss = true;
          	p.setPosition((y_pos + z * Constants::PROTON_MASS_U - nh3_mass)/(DoubleReal)z);
          	p.setIntensity(0.1 / (DoubleReal)(z*z));

						if (z == 1 /*|| y_pos > MIN_DOUBLE_MZ*/)
						{
							spec.push_back(p);
						}
        	}
				}
      }
    }

    // if Q1 abundant loss of water -> pyroglutamic acid formation
		
    if (sequence[0] == 'Q' && prefix == 0 && suffix == 0)
    {
			/*
			for (PeakSpectrum::Iterator it = spec.begin(); it != spec.end(); ++it)
			{
				it->setIntensity(it->getIntensity() * 0.5);
			}*/

			/*		
			for (Size j = 0; j != max_isotope; ++j)
			{
      	p.setPosition((precursor_weight + charge - 1 + j)/(DoubleReal)charge);
      	p.setIntensity(isotope_distributions_[(Int)p.getPosition()[0]][j] * 0.1);
      	spec.push_back(p);
			}
			*/
    }
		

    spec.sortByPosition();

    return;
  }

	void CompNovoIdentificationBase::filterPermuts_(set<String>& permut)
	{
  	set<String> tmp;
  	for (set<String>::const_iterator it = permut.begin(); it != permut.end(); ++it)
  	{
    	if (tryptic_only_)
			{
    		if ((*it)[it->size() - 1] == 'K' || (*it)[it->size() - 1] == 'R')
    		{
      		tmp.insert(*it);
    		}
			}
			else
			{
	    	tmp.insert(*it);
			}
  	}
  	permut = tmp;
  	return;
	}

	Size CompNovoIdentificationBase::countMissedCleavagesTryptic_(const String& peptide) const
	{
  	Size missed_cleavages(0);

  	if (peptide.size() < 2)
  	{
    	return 0;
  	}
  	for (Size i = 0; i != peptide.size() - 1; ++i)
  	{
    	if ((peptide[i] == 'R' || peptide[i] == 'K') && peptide[i + 1] != 'P')
    	{
      	++missed_cleavages;
    	}
  	}

  	return missed_cleavages;
	}

	void CompNovoIdentificationBase::permute_(String prefix, String s, set<String>& permutations)
	{
  	if (s.size() <= 1)
  	{
    	permutations.insert(prefix + s);
  	}
  	else
 		{
    	for (String::Iterator p = s.begin(); p <s.end(); p++ )
    	{
      	char c = *p;
      	s.erase(p);
      	permute_(prefix + c, s , permutations);
      	s.insert(p, c);
    	}
  	}
	}

	void CompNovoIdentificationBase::getDecompositions_(vector<MassDecomposition>& decomps, DoubleReal mass, bool no_caching)
	{
		//static Map<DoubleReal, vector<MassDecomposition> > decomp_cache;
		if (!no_caching)
		{
			if (decomp_cache_.has(mass))
			{
				decomps = decomp_cache_[mass];
				return;
			}
		}

		mass_decomp_algorithm_.getDecompositions(decomps, mass);
		filterDecomps_(decomps);

		if (!no_caching)
		{
			decomp_cache_[mass]  = decomps;
		}
					
  	return;
	}

	void CompNovoIdentificationBase::selectPivotIons_(vector<Size>& pivots, Size left, Size right, Map<DoubleReal, CompNovoIonScoringBase::IonScore>& ion_scores, const PeakSpectrum& CID_spec, DoubleReal precursor_weight, bool full_range)
	{
#ifdef SELECT_PIVOT_DEBUG
  	cerr << "void selectPivotIons(pivots[" << pivots.size() << "], " << left << "[" << CID_spec[left].getPosition()[0] << "]" << ", " << right << "[" << CID_spec[right].getPosition()[0]  << "])" << endl;
#endif
		
		Size max_number_pivot(param_.getValue("max_number_pivot"));

  	// TODO better heuristic, MAX_PIVOT dynamic from range
  	if (right - left > 1)
  	{
    	right -= 1;
    	left += 1;
    	if (right - left < 1 || CID_spec[right].getPosition()[0] - CID_spec[left].getPosition()[0] < 57.0 - fragment_mass_tolerance_)
			{
				return;
			}
			// use more narrow window
			// diff between border and new pivot should be at least 57 - fragment_mass_tolerance (smallest aa)

    	Size new_right(right), new_left(left);
    	for (Size i = left - 1; i < right && CID_spec[i].getPosition()[0] - CID_spec[left - 1].getPosition()[0] < 57.0 - fragment_mass_tolerance_; ++i)
    	{
      	new_left = i;
    	}

    	for (Size i = right + 1; i > new_left &&
         CID_spec[right + 1].getPosition()[0] - CID_spec[i].getPosition()[0] < 57.0 - fragment_mass_tolerance_;
         --i)
    	{
      	new_right = i;
    	}
#ifdef SELECT_PIVOT_DEBUG
			cerr << "new_left=" << new_left << "(" << CID_spec[new_left].getPosition()[0] << "), new_right=" << new_right << "(" << CID_spec[new_right].getPosition()[0]<< ")" << endl;
#endif
    	left = new_left;
    	right = new_right;


    	if (!(right - left > 1))
    	{
      	return;
    	}


    	Size old_num_used(0);
    	set<Size> used_pos;
    	for (Size p = 0; p != min(right - left - 1, max_number_pivot); ++p)
    	{
      	DoubleReal max(0);
      	Size max_pos(0);

      	bool found_pivot(false);
      	for (Size i = left + 1; i < right; ++i)
      	{
					DoubleReal score = ion_scores[CID_spec[i].getPosition()[0]].score;
					DoubleReal position = CID_spec[i].getPosition()[0];
#ifdef SELECT_PIVOT_DEBUG
					cerr << position << " " << precursor_weight << " " << full_range << " " << score;
#endif
					if (score >= max && used_pos.find(i) == used_pos.end())
          {
            // now check if a very similar ion is already selected +/- 3Da
            //bool has_similar(false);
						/*
            for (set<Size>::const_iterator it = used_pos.begin(); it != used_pos.end(); ++it)
            {
             	if (fabs(CID_spec[*it].getPosition()[0] - CID_spec[i].getPosition()[0]) < 1.5)
             	{
               	has_similar = true;
             	}
            }*/

						// TODO this rule should be toggable
						if (!(full_range && (position < precursor_weight / 4.0 || position > precursor_weight / 4.0 * 3.0)))
						{
#ifdef SELECT_PIVOT_DEBUG
							cerr << " max score greater";
#endif
							max = score;
            	max_pos = i;
            	found_pivot = true;
						}
          }
#ifdef SELECT_PIVOT_DEBUG
					cerr << endl;
#endif
      	}

      	used_pos.insert(max_pos);

     		// no pivot ion was added
      	if (!found_pivot || (old_num_used == used_pos.size() && old_num_used != 0))
      	{
        	return;
      	}
      	else
      	{
        	old_num_used = used_pos.size();
      	}

      	pivots.push_back(max_pos);
      	max = 0;
    	}
 		}
  	return;
	}


	// s1 should be the original spectrum
	DoubleReal CompNovoIdentificationBase::compareSpectra_(const PeakSpectrum& s1, const PeakSpectrum& s2)
	{
		DoubleReal score(0.0);
		
		PeakSpectrum::ConstIterator it1 = s1.begin();
		PeakSpectrum::ConstIterator it2 = s2.begin();

		Size num_matches(0);
		while (it1 != s1.end() && it2 != s2.end())
		{
			DoubleReal pos1(it1->getPosition()[0]), pos2(it2->getPosition()[0]);
			if (fabs(pos1 - pos2) < fragment_mass_tolerance_)
			{
				score += it1->getIntensity();
				++num_matches;
			}

			if (pos1 <= pos2)
			{
				++it1;
			}
			else
			{
				++it2;
			}
		}
		
		if (num_matches == 0)
		{
			return 0;
		}
		
		score /= sqrt((DoubleReal)num_matches);
		
		return score;
	}

	
	bool Internal::PermutScoreComparator(const CompNovoIdentificationBase::Permut& p1, const CompNovoIdentificationBase::Permut& p2)
	{
		return p1.getScore() > p2.getScore();
	}

	void CompNovoIdentificationBase::windowMower_(PeakSpectrum& spec, DoubleReal windowsize, Size no_peaks)
	{
  	PeakSpectrum copy(spec);
  	vector<Peak1D> to_be_deleted;
  	for (Size i = 0; i < spec.size(); ++i)
  	{
    	PeakSpectrum sub_spec;
    	bool end(false);
    	for (Size j = i;  spec[j].getPosition()[0] - spec[i].getPosition()[0] < windowsize; )
    	{
      	sub_spec.push_back(spec[j]);
      	if (++j == spec.size())
      	{
        	end = true;
      	  break;
      	}
    	}

    	sub_spec.sortByIntensity(true);

    	for (Size k = no_peaks; k < sub_spec.size(); ++k)
    	{
      	Peak1D p(sub_spec[k]);
      	to_be_deleted.push_back(p);
    	}

    	if (end)
    	{
    	  break;
    	}
  	}

  	spec.clear(false);
  	for (PeakSpectrum::ConstIterator it = copy.begin(); it != copy.end(); ++it)
  	{
    	if (find(to_be_deleted.begin(), to_be_deleted.end(), *it) == to_be_deleted.end())
    	{
      	spec.push_back(*it);
    	}
  	}

  	spec.sortByPosition();

	}


	void CompNovoIdentificationBase::filterDecomps_(vector<MassDecomposition>& decomps)
	{
		Size max_number_aa_per_decomp(param_.getValue("max_number_aa_per_decomp"));
  	vector<MassDecomposition> tmp;
  	for (vector<MassDecomposition>::const_iterator it = decomps.begin(); it != decomps.end(); ++it)
  	{
    	if (it->getNumberOfMaxAA() <= max_number_aa_per_decomp)
    	{
      	tmp.push_back(*it);
    	}
  	}
  	decomps = tmp;
  	return;
	}

	void CompNovoIdentificationBase::initIsotopeDistributions_()
	{
  	IsotopeDistribution iso_dist(max_isotope_);
  	for (Size i = 1; i <= max_mz_ * 2; ++i)
  	{
    	iso_dist.estimateFromPeptideWeight((DoubleReal)i);
    	iso_dist.renormalize();
    	vector<DoubleReal> iso(max_isotope_, 0.0);

    	for (Size j = 0; j != iso_dist.size(); ++j)
    	{
      	iso[j] = iso_dist.getContainer()[j].second;
    	}
    	isotope_distributions_[i] = iso;
  	}
	}

	AASequence CompNovoIdentificationBase::getModifiedAASequence_(const String& sequence)
	{
		AASequence seq;
		for (String::ConstIterator it = sequence.begin(); it != sequence.end(); ++it)
		{
			if (name_to_residue_.has(*it))
			{
				seq += name_to_residue_[*it];
			}
			else
			{
				seq += *it;
			}
		}

		return seq;
	}

	String CompNovoIdentificationBase::getModifiedStringFromAASequence_(const AASequence& sequence)
	{
		String seq;
		for (AASequence::ConstIterator it = sequence.begin(); it != sequence.end(); ++it)
		{
			if (residue_to_name_.has(&*it))
			{
				seq += residue_to_name_[&*it];
			}
			else
			{
				seq += it->getOneLetterCode();
			}
		}
		return seq;
	}

	void 	CompNovoIdentificationBase::updateMembers_()
	{
		// init residue mass table
		String residue_set(param_.getValue("residue_set"));

		set<const Residue*> residues = ResidueDB::getInstance()->getResidues(residue_set);
		for (set<const Residue*>::const_iterator it = residues.begin(); it != residues.end(); ++it)
		{
			aa_to_weight_[(*it)->getOneLetterCode()[0]] = (*it)->getMonoWeight(Residue::Internal);
		}

		max_number_aa_per_decomp_ = param_.getValue("max_number_aa_per_decomp");
		tryptic_only_ = param_.getValue("tryptic_only").toBool();
		fragment_mass_tolerance_ = (DoubleReal)param_.getValue("fragment_mass_tolerance");
		max_number_pivot_ = param_.getValue("max_number_pivot");
		decomp_weights_precision_ = (DoubleReal)param_.getValue("decomp_weights_precision");
		min_mz_ = (DoubleReal)param_.getValue("min_mz");
		max_mz_ = (DoubleReal)param_.getValue("max_mz");
		max_decomp_weight_ = (DoubleReal)param_.getValue("max_decomp_weight");
		max_subscore_number_ = param_.getValue("max_subscore_number");
		max_isotope_ = param_.getValue("max_isotope");

		name_to_residue_.clear();
		residue_to_name_.clear();

		// now handle the modifications
		ModificationDefinitionsSet mod_set((StringList)param_.getValue("fixed_modifications"), (StringList)param_.getValue("variable_modifications"));
		set<ModificationDefinition> fixed_mods = mod_set.getFixedModifications();

		for (set<ModificationDefinition>::const_iterator it = fixed_mods.begin(); it != fixed_mods.end(); ++it)
		{
			ResidueModification mod = ModificationsDB::getInstance()->getModification(it->getModification());
			char aa=' ';
			if (mod.getOrigin().size() != 1 || mod.getOrigin() == "X")
			{
				cerr << "Warning: cannot handle modification " << it->getModification() << ", because aa is ambiguous (" << mod.getOrigin() << "), ignoring modification!" << endl;
				continue;
			}
			else
			{
				aa = mod.getOrigin()[0];
			}

			if (mod.getMonoMass() != 0)
			{
				aa_to_weight_[aa] = mod.getMonoMass();
			}
			else
			{
				if (mod.getDiffMonoMass() != 0)
				{
					aa_to_weight_[aa] += mod.getDiffMonoMass();
				}
				else
				{
					cerr << "Warning: cannot handle modification " << it->getModification() << ", because no monoisotopic mass value was found! Ignoring modification!" << endl;
					continue;
				}
			}

			//cerr << "Setting fixed modification " << it->getModification() << " of amino acid '" << aa << "'; weight = " << aa_to_weight_[aa] << endl;

			const Residue* res = ResidueDB::getInstance()->getModifiedResidue(it->getModification());
			name_to_residue_[aa] = res;
			residue_to_name_[res] = aa;
		}
	
		const StringList mod_names(StringList::create("a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z"));
		vector<String>::const_iterator actual_mod_name = mod_names.begin();
		set<ModificationDefinition> var_mods = mod_set.getVariableModifications();
		for (set<ModificationDefinition>::const_iterator it = var_mods.begin(); it != var_mods.end(); ++it)
		{
			ResidueModification mod = ModificationsDB::getInstance()->getModification(it->getModification());
			char aa = (*actual_mod_name)[0];
			char origin_aa = ' ';
			++actual_mod_name;

			if (mod.getOrigin().size() != 1 || mod.getOrigin() == "X")
			{
				cerr << "CompNovoIdentificationBase: Warning: cannot handle modification " << it->getModification() << ", because aa is ambiguous (" << mod.getOrigin() << "), ignoring modification!" << endl;
				continue;
			}
			else
			{
				origin_aa = mod.getOrigin()[0];
			}

			if (mod.getMonoMass() != 0)
			{
				aa_to_weight_[aa] = mod.getMonoMass();
			}
			else
			{
				if (mod.getDiffMonoMass() != 0)
				{
					aa_to_weight_[aa] = aa_to_weight_[origin_aa] + mod.getDiffMonoMass();
				}
				else
				{
					cerr << "CompNovoIdentificationBase: Warning: cannot handle modification " << it->getModification() << ", because no monoisotopic mass value was found! Ignoring modification!" << endl;
					continue;
				}
			}

			//cerr << "Mapping variable modification " << it->getModification() << " to letter '" << aa << "' (@" << origin_aa << "); weight = " << aa_to_weight_[aa] << endl;
			const Residue* res = ResidueDB::getInstance()->getModifiedResidue(it->getModification());
			name_to_residue_[aa] = res;
			residue_to_name_[res] = aa;
		}

		/*
		cerr << "Following masses are used for identification: " << endl;
		
		for (Map<char, DoubleReal>::const_iterator it = aa_to_weight_.begin(); it != aa_to_weight_.end(); ++it)
		{
			cerr << it->first << " " << precisionWrapper(it->second) << endl;
		}*/
		
		initIsotopeDistributions_();

		Param decomp_param(mass_decomp_algorithm_.getParameters());
		decomp_param.setValue("tolerance", fragment_mass_tolerance_);
		decomp_param.setValue("fixed_modifications", (StringList)param_.getValue("fixed_modifications"));
		decomp_param.setValue("variable_modifications", (StringList)param_.getValue("variable_modifications"));
		mass_decomp_algorithm_.setParameters(decomp_param);
		
		min_aa_weight_ = numeric_limits<DoubleReal>::max();
		for (Map<char, DoubleReal>::const_iterator it = aa_to_weight_.begin(); it != aa_to_weight_.end(); ++it)
		{
			if (min_aa_weight_ > it->second)
			{
				min_aa_weight_ = it->second;
			}
		}
		return;
	}
}

