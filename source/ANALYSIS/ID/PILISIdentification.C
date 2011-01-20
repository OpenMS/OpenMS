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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/ID/PILISIdentification.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/ANALYSIS/ID/PILISScoring.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/KERNEL/StandardTypes.h>

using namespace std;

namespace OpenMS
{
	
	PILISIdentification::PILISIdentification()
		:	DefaultParamHandler("PILISIdentification"),
			/*sequence_db_(0),*/
			hmm_model_(0),
			pre_scorer_(0),
			scorer_(0),
			own_sequence_db_(false),
			own_model_(false)
	{
		defaults_.setValue("precursor_mass_tolerance", 3.0, "Precursor mass tolerance which is used to query the peptide database for peptides");
		defaults_.setValue("peak_mass_tolerance", 0.3, "Peak mass tolerance to align the simulated and experimental spectra");
		defaults_.setValue("max_candidates", 200, "Number of candidates which are kept at the end of the identification");
		defaults_.setValue("pre_score_name", "ZhangSimilarityScore", "The prescoring which is used", StringList::create("advanced"));
		defaults_.setValue("score_name", "ZhangSimilarityScore", "The scoring for the comparison of simulated and experimental spectrum", StringList::create("advanced"));
		defaults_.setValue("use_evalue_scoring", 1, "If set to 1 EValue scoring as described in PILISScoring is used, otherwise similarity scores are directly reported");
		defaults_.setValue("fixed_modifications", "", "fixed modifications to used in the format 57.001@C");
		
		defaultsToParam_();
		updateMembers_();
	}

	PILISIdentification::PILISIdentification(const PILISIdentification& rhs)
		: DefaultParamHandler(rhs),
	/*		sequence_db_(0),*/
			hmm_model_(0),
			own_sequence_db_(false),
			own_model_(false)
	{
		updateMembers_();
	}

	PILISIdentification& PILISIdentification::operator = (const PILISIdentification& rhs)
	{
		if (this != &rhs)
		{
			DefaultParamHandler::operator=(rhs);
		/*	sequence_db_ = 0;*/
			hmm_model_ = 0;
			own_sequence_db_ = false;
			own_model_ = false;
			updateMembers_();
		}
		return *this;
	}

	PILISIdentification::~PILISIdentification()
	{
					/*
		if (own_sequence_db_)
		{
			delete sequence_db_;
		}*/
		if (own_model_)
		{
			delete hmm_model_;
		}
	}
/*
	void PILISIdentification::setSequenceDB(PILISSequenceDB* sequence_db)
	{
		if (own_sequence_db_)
		{
			delete sequence_db_;
			own_sequence_db_ = false;
		}
		sequence_db_ = sequence_db;
	}
*/
	void PILISIdentification::setModel(PILISModel* hmm_model)
	{
		if (own_model_)
		{
			delete hmm_model_;
    	own_model_ = false;
		}
		hmm_model_ = hmm_model;
		return;
	}

	PILISModel* PILISIdentification::getPILISModel_()
	{
		if (hmm_model_ == 0)
		{
			hmm_model_ = new PILISModel();
			own_model_ = true;
		}
		return hmm_model_;
	}
	
	void PILISIdentification::getIdentifications(const vector<map<String, UInt> >& candidates, vector<PeptideIdentification>& ids, const RichPeakMap& exp)
	{
		UInt max_candidates = (UInt)param_.getValue("max_candidates");
		UInt count(0);
		for (RichPeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it, ++count)
		{
			if (it->getMSLevel() != 2)
			{
				continue;
			}

			//cerr << count << "/" << exp.size() << endl;
			PeptideIdentification id;
			getIdentification(candidates[count], id, *it);
			
			//if (id.getHits().size() > max_candidates)
			//{
			//	id.getHits().resize(max_candidates);
			//}

			ids.push_back(id);
		}

		if ((Size)param_.getValue("use_evalue_scoring") != 0)
		{
			PILISScoring scoring;
			scoring.getScores(ids);
		}

		for (vector<PeptideIdentification>::iterator it = ids.begin(); it != ids.end(); ++it)
		{
			if (it->getHits().size() > max_candidates)
			{
				vector<PeptideHit> tmp_hits = it->getHits();
				tmp_hits.resize(max_candidates);
				it->setHits(tmp_hits);
			}
		}

		return;
	}

	void PILISIdentification::getIdentification(const map<String, UInt>& candidates, PeptideIdentification& id, const RichPeakSpectrum& spec)
	{
		if (spec.getMSLevel() != 2)
		{
			return;
		}
		
		RichPeakSpectrum spec_copy(spec);
		Normalizer normalizer;
		Param param(normalizer.getParameters());
		param.setValue("method", "to_one");
		normalizer.setParameters(param);

		normalizer.filterSpectrum(spec_copy);

		//double pre_tol = (double)param_.getValue("precursor_mass_tolerance");
		String score_name = param_.getValue("score_name");
		
		scorer_ = Factory<PeakSpectrumCompareFunctor>::create(score_name);
		Param scorer_param(scorer_->getParameters());
		scorer_param.setValue("epsilon",(DoubleReal)param_.getValue("peak_mass_tolerance"));
		scorer_->setParameters(scorer_param);
	
		double pre_pos = 0.0;
		if (!spec_copy.getPrecursors().empty()) pre_pos = spec_copy.getPrecursors()[0].getMZ();
    if (pre_pos < 200) // TODO
    {
      cerr << "PILISIdentification: spectrum does not have a precursor peak set. Precursor peak @ m/z=" << pre_pos << endl;
      return;
    }


		//cerr << "#cand peptides: " << cand_peptides.size() << ", " << pre_pos << ", +/- " << pre_tol << endl;

		PeptideIdentification pre_id;
		getPreIdentification_(pre_id, spec_copy, candidates);

		getFinalIdentification_(id, spec_copy, pre_id);
/*
		SpectrumAlignment aligner;
		Param aligner_param(aligner.getParameters());
		aligner_param.setValue("epsilon", 0.3);
		aligner.setParameters(aligner_param);

		for (Size i = 0; i != id.getPeptideHits().size(); ++i)
		{
			vector<pair<UInt, UInt> > alignment;
			aligner.getSpectrumAlignment(alignment, spec_copy, sim_specs_[i]);

			cerr << i << " " << id.getPeptideHits()[i].getSequence() << endl;
			double rms(0);
			for (Size j = 0; j != alignment.size(); ++j)
			{
				double mz1(spec_copy[alignment[j].first].getMZ()), mz2(sim_specs_[i][alignment[j].second].getMZ());
				cerr << mz1 << " " << mz2 << " " << mz1 - mz2 << " " << sim_specs_[i][alignment[j].second].getMetaValue("IonName") << endl;
				rms += pow(mz1 - mz2, 2.0);
			}
			cerr << "RMS=" << sqrt(rms/double(id.getPeptideHits().size())) << endl;
		}*/

		if ((Size)param_.getValue("use_evalue_scoring") != 0)
		{
			PILISScoring scoring;
			scoring.getScore(id);
		}

		UInt max_candidates = (UInt)param_.getValue("max_candidates");
		if (id.getHits().size() > max_candidates)
		{
			vector<PeptideHit> tmp_hits = id.getHits();
			tmp_hits.resize(max_candidates);
			id.setHits(tmp_hits);
		}

		

		return;
	}

	void PILISIdentification::getPreIdentification_(PeptideIdentification& id, const RichPeakSpectrum& spec, const map<String, UInt>& cand_peptides)
	{
    // get simple spectra to pre-eliminate most of the candidates
    for (map<String, UInt>::const_iterator it1 = cand_peptides.begin(); it1 != cand_peptides.end(); ++it1)
    {
      // TODO parameter settings
      RichPeakSpectrum sim_spec;
			//cerr << it1->first << " " << it1->second << endl;
			try 
			{
				AASequence seq(it1->first);
      	getSpectrum_(sim_spec, it1->first, it1->second);
			}
			catch (Exception::ParseError e)
			{
				cerr << "Peptide sequence " << it1->first << " cannot be processed" << endl;
				continue;
			}

			//TODO WARNING ERROR KOTZ (Andreas, Marc)
			PeakSpectrum s1,s2;
			s1.resize(sim_spec.size());
			for (Size p=0; p<sim_spec.size(); ++p)
			{
				s1[p] = sim_spec[p];
			}
			s2.resize(spec.size());
			for (Size p=0; p< spec.size(); ++p)
			{
				s2[p] = spec[p];
			}
      double score = (*scorer_)(s1, s2);
			//cerr << "Pre: " << it1->first << " " << it1->second << " " << score << endl;
      PeptideHit peptide_hit(score, 0, it1->second, it1->first);
      id.insertHit(peptide_hit);
    }

    id.assignRanks();
		return;
	}

	void PILISIdentification::getFinalIdentification_(PeptideIdentification& id, const RichPeakSpectrum& spec, const PeptideIdentification& pre_id)
	{
		UInt max_candidates = (UInt)param_.getValue("max_candidates");
		sim_specs_.clear();
		id.setScoreType("PILIS");
		for (Size i = 0; i < pre_id.getHits().size() && i < max_candidates; ++i)
    {
      AASequence peptide_sequence = pre_id.getHits()[i].getSequence();
      RichPeakSpectrum sim_spec;
      getPILISModel_()->getSpectrum(sim_spec, peptide_sequence, pre_id.getHits()[i].getCharge());
			sim_specs_.push_back(sim_spec);
			
			//TODO WARNING ERROR KOTZ (Andreas, Marc)
			PeakSpectrum s1,s2;
			s1.resize(sim_spec.size());
			for (Size p=0; p<sim_spec.size(); ++p)
			{
				s1[p] = sim_spec[p];
			}
			s2.resize(spec.size());
			for (Size p=0; p< spec.size(); ++p)
			{
				s2[p] = spec[p];
			}
      double score = (*scorer_)(s1, s2);
			//cerr << "Final: " << peptide_sequence << " " << pre_id.getHits()[i].getCharge() << " " << score << endl;
      PeptideHit peptide_hit(score, 0, pre_id.getHits()[i].getCharge(), peptide_sequence);
      id.insertHit(peptide_hit);
    }

    id.assignRanks();

		return;
	}
	
	void PILISIdentification::getSpectrum_(RichPeakSpectrum& spec, const String& sequence, int charge)
	{
		double b_pos(0);
		double y_pos(18);
		bool b_H2O_loss(false), b_NH3_loss(false), y_H2O_loss(false), y_NH3_loss(false);
		for (Size i = 0; i != sequence.size(); ++i)
		{
			char aa(sequence[i]);
			b_pos += aa_weight_[aa];
			
			char aa2(sequence[sequence.size() - i - 1]);
			y_pos += aa_weight_[aa2];
			for (int z = 1; z <= charge && z < 3; ++z)
			{
				// b-ions
				p_.setPosition((b_pos + z)/z);
				p_.setIntensity(0.8f);
				spec.push_back(p_);

				// b-ion losses
				if (b_H2O_loss || aa == 'S' || aa == 'T' || aa == 'E' || aa == 'D')
				{
					b_H2O_loss = true;
					p_.setPosition((b_pos + z - 18.0)/z);
					p_.setIntensity(0.1f);
					spec.push_back(p_);
				}
				if (b_NH3_loss || aa == 'Q' || aa == 'N' || aa == 'R' || aa == 'K')
				{
					b_NH3_loss = true;
					p_.setPosition((b_pos + z - 17.0)/z);
					p_.setIntensity(0.1f);
					spec.push_back(p_);
				}

				// a-ions
				p_.setPosition((b_pos +z - 28.0)/z);
				p_.setIntensity(0.3f);
				spec.push_back(p_);
				
				// y-ions
				p_.setPosition((y_pos + z)/z);
				p_.setIntensity(1.0);
				spec.push_back(p_);

				if (y_H2O_loss || aa2 == 'S' || aa2 == 'T' || aa2 == 'E' || aa2 == 'D' || aa2 == 'Q')
				{
					y_H2O_loss = true;
					p_.setPosition((y_pos + z - 18.0)/z);
					if (aa2 != 'Q')
					{
						p_.setIntensity(0.2f);
					}
					else
					{
						p_.setIntensity(1);
					}
					spec.push_back(p_);
				}
				if (y_NH3_loss || aa == 'Q' || aa == 'N' || aa == 'R' || aa == 'K')
				{
					y_NH3_loss = true;
					p_.setPosition((y_pos + z - 17.0)/z);
					p_.setIntensity(0.2f);
				}
			}
		}

		spec.sortByPosition();
		
		return;
	}

	void PILISIdentification::updateMembers_()
	{
		pre_scorer_ = Factory<PeakSpectrumCompareFunctor>::create((String)defaults_.getValue("pre_score_name"));
    scorer_ = Factory<PeakSpectrumCompareFunctor>::create((String)defaults_.getValue("score_name"));

		// set amino acids to the default weights
    aa_weight_['K'] = 128.095;
    aa_weight_['M'] = 131.04;
    aa_weight_['F'] = 147.068;
    aa_weight_['P'] = 97.0528;
    aa_weight_['S'] = 87.032;
    aa_weight_['T'] = 101.048;
    aa_weight_['W'] = 186.079;
    aa_weight_['Y'] = 163.063;
    aa_weight_['V'] = 99.0684;
    aa_weight_['A'] = 71.0371;
    aa_weight_['R'] = 156.101;
    aa_weight_['N'] = 114.043;
    aa_weight_['D'] = 115.027;
    //aa_weight_['C'] = 161.015; //CmC
    aa_weight_['C'] = 103.00919;
    aa_weight_['E'] = 129.043;
    aa_weight_['Q'] = 128.059;
    aa_weight_['G'] = 57.0215;
    aa_weight_['H'] = 137.059;
    aa_weight_['I'] = 113.084;
    aa_weight_['L'] = 113.084;

		// decode the fixed modifications string
		String fixed_modifications = param_.getValue("fixed_modifications");

		if (fixed_modifications != "")
		{
			//cerr << fixed_modifications << endl;
			vector<String> mod_split;

			fixed_modifications.split(',', mod_split); // comma separated modifications

			// now get the modifications
			for (Size i = 0; i != mod_split.size(); ++i)
			{
				for (Size j = 0; j != mod_split[i].size(); ++j)
				{
					if (mod_split[i][j] == '@')
					{
						DoubleReal mass_diff(mod_split[i].substr(0, j).toFloat());
						if (j != mod_split[i].size() - 2)
						{
							throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "modification in wrong format", "weight@one_letter_code");
						}
						char res = mod_split[i][j + 1];

						if (aa_weight_.has(res))
						{
							//cerr << res << " " << mass_diff << endl;
							aa_weight_[res] += mass_diff;
						}
					}
				}
			}
		}
	}
}

