// -*- Mode: C++; tab-width: 2; -*-
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
//

#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationCID.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FILTERING/TRANSFORMERS/Normalizer.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringCID.h>

#include <boost/math/special_functions/fpclassify.hpp>

//#define DAC_DEBUG

//#define WRITE_SCORED_SPEC
//#define REDUCE_PERMUTS_DEBUG

#ifdef WRITE_SCORED_SPEC
	#include <OpenMS/FORMAT/DTAFile.h>
#endif

//#define SPIKE_IN

using namespace std;

namespace OpenMS
{
	CompNovoIdentificationCID::CompNovoIdentificationCID()
		:	CompNovoIdentificationBase(),
			precursor_mass_tolerance_(1.5)
	{
		defaultsToParam_();
		updateMembers_();
	}

	CompNovoIdentificationCID::CompNovoIdentificationCID(const CompNovoIdentificationCID& rhs)
		: CompNovoIdentificationBase(rhs)
	{
		updateMembers_();
	}

	CompNovoIdentificationCID& CompNovoIdentificationCID::operator = (const CompNovoIdentificationCID& rhs)
	{
		if (this != &rhs)
		{
			CompNovoIdentificationBase::operator = (rhs);
			updateMembers_();
		}
		return *this;
	}

	CompNovoIdentificationCID::~CompNovoIdentificationCID()
	{
	}

	void CompNovoIdentificationCID::getIdentifications(vector<PeptideIdentification>& pep_ids, const PeakMap& exp)
	{
		UInt count(1);
		for (PeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it, ++count)
		{
			cerr << count << "/" << exp.size() << endl;
			PeptideIdentification id;
			// TODO check if both CID and ETD is present;
			PeakSpectrum CID_spec(*it);
			id.setMetaValue("RT", it->getRT());
			id.setMetaValue("MZ", it->getPrecursors().begin()->getMZ());
			
			subspec_to_sequences_.clear();
			permute_cache_.clear();
			decomp_cache_.clear();
			
			getIdentification(id, CID_spec);
			cerr << "size_of id=" << id.getHits().size() << endl;
			pep_ids.push_back(id);

			//++it;
			
			//
			//if (count == 10)
			//{
			//	return;
			//}
		}
		return;
	}

	void CompNovoIdentificationCID::getIdentification(PeptideIdentification& id, const PeakSpectrum& CID_spec)
	{
		//if (CID_spec.getPrecursors().begin()->getMZ() > 1000.0)
		//{
		//	cerr << "Weight of precursor has been estimated to exceed 2000.0 Da which is the current limit" << endl;
		//	return;
		//}
					
		PeakSpectrum new_CID_spec(CID_spec);
		windowMower_(new_CID_spec, 0.3, 1);

		Param zhang_param;
		zhang_param = zhang_.getParameters();
		zhang_param.setValue("tolerance", fragment_mass_tolerance_);
		zhang_param.setValue("use_linear_factor", "true");
		zhang_.setParameters(zhang_param);
								
		
		Normalizer normalizer;
  	Param n_param(normalizer.getParameters());
  	n_param.setValue("method", "to_one");
  	normalizer.setParameters(n_param);
  	normalizer.filterSpectrum(new_CID_spec);

		UInt charge(2);
		DoubleReal precursor_weight(0); // [M+H]+
		if (CID_spec.getPrecursors().size() > 0)
		{
			// believe charge of spectrum?
			if (CID_spec.getPrecursors().begin()->getCharge() != 0)
			{
				charge = CID_spec.getPrecursors().begin()->getCharge();
			}
			else
			{
				// TODO estimate charge state
			}
			precursor_weight = CID_spec.getPrecursors().begin()->getMZ() * charge - ((charge - 1) * Constants::PROTON_MASS_U);
		}

		cerr << "charge=" << charge << ", [M+H]=" << precursor_weight << endl;
				
		// now delete all peaks that are right of the estimated precursor weight
		UInt peak_counter(0);
		for (PeakSpectrum::ConstIterator it = new_CID_spec.begin(); it != new_CID_spec.end(); ++it, ++peak_counter)
		{
			if (it->getPosition()[0] > precursor_weight)
			{
				break;
			}
		}
		if (peak_counter < new_CID_spec.size())
		{
			new_CID_spec.resize(peak_counter);
		}

		
		static DoubleReal oxonium_mass = EmpiricalFormula("H2O+").getMonoWeight();

		Peak1D p;
		p.setIntensity(1);
		p.setPosition(oxonium_mass);

		new_CID_spec.push_back(p);
		
		p.setPosition(precursor_weight);
		new_CID_spec.push_back(p);
		
		// add complement to spectrum
		/*
    for (PeakSpectrum::ConstIterator it1 = CID_spec.begin(); it1 != CID_spec.end(); ++it1)
    {
      // get m/z of complement
      DoubleReal mz_comp = precursor_weight - it1->getPosition()[0] + Constants::PROTON_MASS_U;

      // search if peaks are available that have similar m/z values
      UInt count(0);
      bool found(false);
      for (PeakSpectrum::ConstIterator it2 = CID_spec.begin(); it2 != CID_spec.end(); ++it2, ++count)
      {
        if (fabs(mz_comp - it2->getPosition()[0]) < fragment_mass_tolerance)
        {
          // add peak intensity to corresponding peak in new_CID_spec
          new_CID_spec[count].setIntensity(new_CID_spec[count].getIntensity());
        }
      }
      if (!found)
      {
        // infere this peak
        Peak1D p;
        p.setIntensity(it1->getIntensity());
        p.setPosition(mz_comp);
        new_CID_spec.push_back(p);
      }
    }*/
	
		CompNovoIonScoringCID ion_scoring;
		Param ion_scoring_param(ion_scoring.getParameters());
		ion_scoring_param.setValue("fragment_mass_tolerance", fragment_mass_tolerance_);
		ion_scoring_param.setValue("precursor_mass_tolerance", precursor_mass_tolerance_);
		ion_scoring_param.setValue("decomp_weights_precision", decomp_weights_precision_);
		ion_scoring_param.setValue("double_charged_iso_threshold", (DoubleReal)param_.getValue("double_charged_iso_threshold"));
		ion_scoring_param.setValue("max_isotope_to_score", (UInt)param_.getValue("max_isotope_to_score"));
		ion_scoring_param.setValue("max_isotope", max_isotope_);
		ion_scoring.setParameters(ion_scoring_param);
		
		Map<DoubleReal, IonScore> ion_scores;
		ion_scoring.scoreSpectrum(ion_scores, new_CID_spec, precursor_weight, charge);

		new_CID_spec.sortByPosition();
		
		/*
		cerr << "Size of ion_scores " << ion_scores.size() << endl;
		for (Map<DoubleReal, IonScore>::const_iterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
		{
			cerr << it->first << " " << it->second.score << endl;
		}*/
		
#ifdef WRITE_SCORED_SPEC
		PeakSpectrum filtered_spec(new_CID_spec);
  	filtered_spec.clear();
  	for (Map<DoubleReal, CompNovoIonScoringCID::IonScore>::const_iterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
  	{
    	Peak1D p;
    	p.setIntensity(it->second.score);
    	p.setPosition(it->first);
    	filtered_spec.push_back(p);
  	}
  	DTAFile().store("spec_scored.dta", filtered_spec);
#endif

		set<String> sequences;
		getDecompositionsDAC_(sequences, 0, new_CID_spec.size() - 1, precursor_weight, new_CID_spec, ion_scores);

#ifdef SPIKE_IN
		sequences.insert("AFCVDGEGR");
		sequences.insert("APEFAAPWPDFVPR");
		sequences.insert("AVKQFEESQGR");
		sequences.insert("CCTESLVNR");
		sequences.insert("DAFLGSFLYEYSR");
		sequences.insert("DAIPENLPPLTADFAEDK");
		sequences.insert("DDNKVEDIWSFLSK");
		sequences.insert("DDPHACYSTVFDK");
		sequences.insert("DEYELLCLDGSR");
		sequences.insert("DGAESYKELSVLLPNR");
    sequences.insert("DGASCWCVDADGR");
    sequences.insert("DLFIPTCLETGEFAR");
    sequences.insert("DTHKSEIAHR");
    sequences.insert("DVCKNYQEAK");
    sequences.insert("EACFAVEGPK");
    sequences.insert("ECCHGDLLECADDR");
    sequences.insert("EFLGDKFYTVISSLK");
    sequences.insert("EFTPVLQADFQK");
    sequences.insert("ELFLDSGIFQPMLQGR");
    sequences.insert("ETYGDMADCCEK");
    sequences.insert("EVGCPSSSVQEMVSCLR");
    sequences.insert("EYEATLEECCAK");
    sequences.insert("FADLIQSGTFQLHLDSK");
    sequences.insert("FFSASCVPGATIEQK");
    sequences.insert("FLANVSTVLTSK");
    sequences.insert("FLSGSDYAIR");
    sequences.insert("FTASCPPSIK");
    sequences.insert("GAIEWEGIESGSVEQAVAK");
    sequences.insert("GDVAFIQHSTVEENTGGK");
    sequences.insert("GEPPSCAEDQSCPSER");
    sequences.insert("GEYVPTSLTAR");
    sequences.insert("GQEFTITGQKR");
    sequences.insert("GTFAALSELHCDK");
    sequences.insert("HLVDEPQNLIK");
    sequences.insert("HQDCLVTTLQTQPGAVR");
    sequences.insert("HTTVNENAPDQK");
    sequences.insert("ILDCGSPDTEVR");
    sequences.insert("KCPSPCQLQAER");
    sequences.insert("KGTEFTVNDLQGK");
    sequences.insert("KQTALVELLK");
    sequences.insert("KVPQVSTPTLVEVSR");
    sequences.insert("LALQFTTNAKR");
    sequences.insert("LCVLHEKTPVSEK");
    sequences.insert("LFTFHADICTLPDTEK");
    sequences.insert("LGEYGFQNALIVR");
    sequences.insert("LHVDPENFK");
    sequences.insert("LKECCDKPLLEK");
    sequences.insert("LKHLVDEPQNLIK");
    sequences.insert("LKPDPNTLCDEFK");
    sequences.insert("LLGNVLVVVLAR");
    sequences.insert("LLVVYPWTQR");
    sequences.insert("LRVDPVNFK");
    sequences.insert("LTDEELAFPPLSPSR");
    sequences.insert("LVNELTEFAK");
    sequences.insert("MFLSFPTTK");
    sequences.insert("MPCTEDYLSLILNR");
    sequences.insert("NAPYSGYSGAFHCLK");
    sequences.insert("NECFLSHKDDSPDLPK");
    sequences.insert("NEPNKVPACPGSCEEVK");
    sequences.insert("NLQMDDFELLCTDGR");
    sequences.insert("QAGVQAEPSPK");
    sequences.insert("RAPEFAAPWPDFVPR");
    sequences.insert("RHPEYAVSVLLR");
    sequences.insert("RPCFSALTPDETYVPK");
    sequences.insert("RSLLLAPEEGPVSQR");
    sequences.insert("SAFPPEPLLCSVQR");
    sequences.insert("SAGWNIPIGTLLHR");
    sequences.insert("SCWCVDEAGQK");
    sequences.insert("SGNPNYPHEFSR");
    sequences.insert("SHCIAEVEK");
    sequences.insert("SISSGFFECER");
    sequences.insert("SKYLASASTMDHAR");
    sequences.insert("SLHTLFGDELCK");
    sequences.insert("SLLLAPEEGPVSQR");
    sequences.insert("SPPQCSPDGAFRPVQCK");
    sequences.insert("SREGDPLAVYLK");
    sequences.insert("SRQIPQCPTSCER");
    sequences.insert("TAGTPVSIPVCDDSSVK");
    sequences.insert("TCVADESHAGCEK");
    sequences.insert("TQFGCLEGFGR");
    sequences.insert("TVMENFVAFVDK");
    sequences.insert("TYFPHFDLSHGSAQVK");
    sequences.insert("TYMLAFDVNDEK");
    sequences.insert("VDEVGGEALGR");
    sequences.insert("VDLLIGSSQDDGLINR");
    sequences.insert("VEDIWSFLSK");
    sequences.insert("VGGHAAEYGAEALER");
    sequences.insert("VGTRCCTKPESER");
    sequences.insert("VKVDEVGGEALGR");
    sequences.insert("VKVDLLIGSSQDDGLINR");
    sequences.insert("VLDSFSNGMK");
    sequences.insert("VLSAADKGNVK");
    sequences.insert("VPQVSTPTLVEVSR");
    sequences.insert("VTKCCTESLVNR");
    sequences.insert("VVAASDASQDALGCVK");
    sequences.insert("VVAGVANALAHR");
    sequences.insert("YICDNQDTISSK");
    sequences.insert("YLASASTMDHAR");
    sequences.insert("YNGVFQECCQAEDK");
#endif
		
		SpectrumAlignmentScore spectra_zhang;
		spectra_zhang.setParameters(zhang_param);
	
    vector<PeptideHit> hits;
		UInt number_missed_cleavages = (UInt)param_.getValue("number_missed_cleavages");
    for (set<String>::const_iterator it = sequences.begin(); it != sequences.end(); ++it)
    {
			
			UInt num_missed = countMissedCleavagesTryptic_(*it);
			if (number_missed_cleavages < num_missed)
			{
				cerr << "Two many missed cleavages: " << *it << ", found " << num_missed << ", allowed " << number_missed_cleavages << endl;
				continue;
			}
      PeakSpectrum CID_sim_spec;			
      getCIDSpectrum_(CID_sim_spec, *it, charge);

      //normalizer.filterSpectrum(CID_sim_spec);

      DoubleReal cid_score = zhang_(CID_sim_spec, CID_spec);

      PeptideHit hit;
      hit.setScore(cid_score);

      hit.setSequence(getModifiedAASequence_(*it));
      hit.setCharge(charge);
      hits.push_back(hit);
      cerr << getModifiedAASequence_(*it) << " " << cid_score << " " << endl;
    }
    
		// rescore the top hits
		id.setHits(hits);
		id.assignRanks();
	
		hits = id.getHits();
	
		SpectrumAlignmentScore alignment_score;
		Param align_param(alignment_score.getParameters());
		align_param.setValue("tolerance", fragment_mass_tolerance_);
		align_param.setValue("use_linear_factor", "true");
		alignment_score.setParameters(align_param);

		for (vector<PeptideHit>::iterator it = hits.begin(); it != hits.end(); ++it)
		{
			cerr << "Pre: " << it->getRank() << " " << it->getSequence() << " " << it->getScore() << " " << endl;
		}
		
		UInt number_of_prescoring_hits = (UInt)param_.getValue("number_of_prescoring_hits");
		if (hits.size() > number_of_prescoring_hits)
		{
			hits.resize(number_of_prescoring_hits);
		}

		for (vector<PeptideHit>::iterator it = hits.begin(); it != hits.end(); ++it)
		{
			PeakSpectrum CID_sim_spec;
			getCIDSpectrum_(CID_sim_spec, getModifiedStringFromAASequence_(it->getSequence()), charge);

			normalizer.filterSpectrum(CID_sim_spec);

			//DTAFile().store("sim_specs/" + it->getSequence().toUnmodifiedString() + "_sim_CID.dta", CID_sim_spec);
			
			//DoubleReal cid_score = spectra_zhang(CID_sim_spec, CID_spec);
			DoubleReal cid_score = alignment_score(CID_sim_spec, CID_spec);

			cerr << "Final: " << it->getSequence() << " " << cid_score << endl;
		
			it->setScore(cid_score);
		}

		id.setHits(hits);
		id.assignRanks();
		hits = id.getHits();
		
		for (vector<PeptideHit>::iterator it = hits.begin(); it != hits.end(); ++it)
		{
			cerr << "Fin: " << it->getRank() << " " << it->getSequence() << " " << it->getScore() << " " << endl;
		}

		UInt number_of_hits = (UInt)param_.getValue("number_of_hits");
		if (id.getHits().size() > number_of_hits)
		{
			hits.resize(number_of_hits);
		}
	
		id.setHits(hits);
		id.assignRanks();

		return;
	}

	void CompNovoIdentificationCID::reducePermuts_(set<String>& permuts, const PeakSpectrum& CID_spec, DoubleReal prefix, DoubleReal suffix)
	{
  	if (permuts.size() < max_subscore_number_)
  	{
   		return;
  	}

		vector<Permut> score_permuts;
		
  	UInt i(0);
  	for (set<String>::const_iterator it = permuts.begin(); it != permuts.end(); ++it, ++i)
  	{
#ifdef REDUCE_PERMUTS_DEBUG
    	if (i % 1000 == 0)
    	{
      	cerr << (DoubleReal)i / permuts.size() * 100 << "%" << endl;
    	}
#endif
			
    	PeakSpectrum CID_sim_spec;
    	getCIDSpectrumLight_(CID_sim_spec, *it, prefix, suffix);
			//getCIDSpectrum_(CID_sim_spec, *it, 1, prefix, suffix);
    	
    	DoubleReal score = zhang_(CID_sim_spec, CID_spec);;

			if (boost::math::isnan(score))
			{
				score = 0;
			}

			score /= it->size();

			if (boost::math::isnan(score))
			{
				score = 0;
			}


#ifdef REDUCE_PERMUTS_DEBUG
    	cerr << "Subscoring: " << *it << " " << cid_score << " (CID=";
/*    	for (PeakSpectrum::ConstIterator pit = CID_sim_spec.begin(); pit != CID_sim_spec.end(); ++pit)
    	{
      	cerr << pit->getPosition()[0] << "|" << pit->getIntensity() << "; ";
    	}*/
    	cerr << endl;
#endif

			Permut new_permut(it, score);
			score_permuts.push_back(new_permut);
  	}

		sort(score_permuts.begin(), score_permuts.end(), Internal::PermutScoreComparator);

  	set<String> new_permuts;
  	UInt count(0);
		for (vector<Permut>::const_iterator it = score_permuts.begin(); it != score_permuts.end() && count < max_subscore_number_; ++it, ++count)
		{
			new_permuts.insert(*it->getPermut());
#ifdef REDUCE_PERMUTS_DEBUG
			cerr << "Subscore winner: " << it->getPermut() << " " << it->getScore() << endl;
#endif
		}

  	permuts = new_permuts;
  	return;
	}


// divide and conquer algorithm of the sequencing
void CompNovoIdentificationCID::getDecompositionsDAC_(set<String>& sequences, UInt left, UInt right, DoubleReal peptide_weight, const PeakSpectrum& CID_spec, Map<DoubleReal, CompNovoIonScoringCID::IonScore>& ion_scores)
{
	static DoubleReal oxonium_mass = EmpiricalFormula("H2O+").getMonoWeight();
	DoubleReal offset_suffix(CID_spec[left].getPosition()[0] - oxonium_mass);
	DoubleReal offset_prefix(peptide_weight - CID_spec[right].getPosition()[0]);

#ifdef DAC_DEBUG
	static Int depth_(0);
	++depth_;
	String tabs_(depth_, '\t');
	cerr << tabs_ << "void getDecompositionsDAC(sequences[" << sequences.size() << "], " << left << ", " << right << ") ";
	cerr << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " diff=";
#endif
	
	DoubleReal diff = CID_spec[right].getPosition()[0] - CID_spec[left].getPosition()[0];

#ifdef DAC_DEBUG
	cerr << diff << endl;
	cerr << "offset_prefix=" << offset_prefix << ", offset_suffix=" << offset_suffix << endl;
#endif

	if (subspec_to_sequences_.has(left) && subspec_to_sequences_[left].has(right))
	{
		sequences = subspec_to_sequences_[left][right];

#ifdef DAC_DEBUG
		depth_--;
		cerr << tabs_ << "from cache DAC: " << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << sequences.size() << " " << left << " " << right << endl;
#endif
		return;
	}
	
	// no further solutions possible?
	if (diff < min_aa_weight_)
  {
#ifdef DAC_DEBUG
		depth_--;
#endif	
    return;
  }

	// no further division needed?
  if (diff <= max_decomp_weight_)
  {
		vector<MassDecomposition> decomps;

		// if we are at the C-terminus use precursor_mass_tolerance_
		if (offset_prefix < precursor_mass_tolerance_)
		{
			Param decomp_param(mass_decomp_algorithm_.getParameters());
			decomp_param.setValue("tolerance", precursor_mass_tolerance_);
			mass_decomp_algorithm_.setParameters(decomp_param);
			getDecompositions_(decomps, diff);
			decomp_param.setValue("tolerance", fragment_mass_tolerance_);
			mass_decomp_algorithm_.setParameters(decomp_param);
		}
		else
		{
			getDecompositions_(decomps, diff);
		}
		//filterDecomps_(decomps);
		
#ifdef DAC_DEBUG
		cerr << tabs_ << "Found " << decomps.size() << " decomps" << endl;
		cerr << tabs_ << "Permuting...";
#endif
		
		//static Map<String, set<String> > permute_cache;
		for (vector<MassDecomposition>::const_iterator it = decomps.begin(); it != decomps.end(); ++it)
		{
#ifdef DAC_DEBUG
			cerr << it->toString() << endl;
#endif
		
			String exp_string = it->toExpandedString();
			if (!permute_cache_.has(exp_string))
			{
				permute_("", exp_string, sequences);
				permute_cache_[exp_string] = sequences;
			}
			else
			{
				sequences = permute_cache_[exp_string];
			}
		}
		
#ifdef DAC_DEBUG
		cerr << tabs_ << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << peptide_weight << endl;
		if (sequences.size() > max_subscore_number_)
		{
    	cerr << tabs_ << "Reducing #sequences from " << sequences.size() << " to " << max_subscore_number_ << "(prefix=" << offset_prefix  << ", suffix=" << offset_suffix << ")...";
		}
#endif

		// C-terminus
		if (offset_suffix <= precursor_mass_tolerance_)
		{
			filterPermuts_(sequences);
		}
	
		// reduce the sequences
    reducePermuts_(sequences, CID_spec, offset_prefix, offset_suffix);
#ifdef DAC_DEBUG
		cerr << "Writing to cache " << left << " " << right << endl;
#endif
		subspec_to_sequences_[left][right] = sequences;

#ifdef DAC_DEBUG		
		cerr << "ended" << endl;
		cerr << tabs_ << "DAC: " << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << sequences.size() << endl;
		depth_--;
#endif
		
    return;
  }

	// select suitable pivot peaks
  vector<UInt> pivots;

	if (offset_suffix < precursor_mass_tolerance_ && offset_prefix < precursor_mass_tolerance_)
	{
  	selectPivotIons_(pivots, left, right, ion_scores, CID_spec, peptide_weight, true);
	}
	else
	{
		selectPivotIons_(pivots, left, right, ion_scores, CID_spec, peptide_weight, false);
	}

	// run divide step
#ifdef DAC_DEBUG
	cerr << tabs_ << "Selected " << pivots.size() << " pivot ions: ";
	for (vector<UInt>::const_iterator it = pivots.begin(); it != pivots.end(); ++it)
	{
		cerr << *it << "(" << CID_spec[*it].getPosition()[0] << ") ";
	}
	cerr << endl;
#endif	

  for (vector<UInt>::const_iterator it = pivots.begin(); it != pivots.end(); ++it)
  {
		set<String> seq1, seq2, new_sequences;
		
		// the smaller the 'gap' the greater the chance of not finding anything
		// so we we compute the smaller gap first
		DoubleReal diff1(CID_spec[*it].getPosition()[0] - CID_spec[left].getPosition()[0]);
		DoubleReal diff2(CID_spec[right].getPosition()[0] - CID_spec[*it].getPosition()[0]);
		
		if (diff1 < diff2)
		{
    	getDecompositionsDAC_(seq1, left, *it, peptide_weight, CID_spec, ion_scores);
			if (seq1.size() == 0)
			{
#ifdef DAC_DEBUG
				cerr << tabs_ << "first call produced 0 candidates (" << diff1 << ")" << endl;
#endif
				continue;
			}

			getDecompositionsDAC_(seq2, *it, right, peptide_weight, CID_spec, ion_scores);
		}
		else
		{
			getDecompositionsDAC_(seq2, *it, right, peptide_weight, CID_spec, ion_scores);
			if (seq2.size() == 0)
			{
#ifdef DAC_DEBUG
				cerr << tabs_ << "second call produced 0 candidates (" << diff2 << ")" << endl;
#endif
				continue;
			}

			getDecompositionsDAC_(seq1, left, *it, peptide_weight, CID_spec, ion_scores);
		}

#ifdef DAC_DEBUG
		cerr << tabs_ << "Found " << seq1.size() << " solutions (1) " << diff1 << endl;
		cerr << tabs_ << "Found " << seq2.size() << " solutions (2) " << diff2 << endl;
		cerr << tabs_ << "inserting " << seq1.size() * seq2.size()  << " sequences" << endl;
#endif

		// C-terminus
		if (offset_suffix <= fragment_mass_tolerance_)
		{
			filterPermuts_(seq1);
		}
		
		// test if we found enough sequence candidates
		if (seq1.size() == 0 || seq2.size() == 0)
		{
			continue;
		}
		
    for (set<String>::const_iterator it1 = seq1.begin(); it1 != seq1.end(); ++it1)
    {
      for (set<String>::const_iterator it2 = seq2.begin(); it2 != seq2.end(); ++it2)
     	{
				new_sequences.insert(*it2 + *it1);
      }
    }

		if (seq1.size() * seq2.size() > max_subscore_number_/* && (offset_prefix > fragment_mass_tolerance_ || offset_suffix > fragment_mass_tolerance_)*/)
		{
#ifdef DAC_DEBUG			
			cerr << tabs_ << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << peptide_weight << endl;
			cerr << tabs_ << "Reducing #sequences from " << new_sequences.size() << " to " << max_subscore_number_ << "(prefix=" << offset_prefix  << ", suffix=" << offset_suffix << ")...";
#endif
			if (offset_prefix > precursor_mass_tolerance_ || offset_prefix > precursor_mass_tolerance_)
			{
				reducePermuts_(new_sequences, CID_spec, offset_prefix, offset_suffix);
			}
			
#ifdef DAC_DEBUG
			for (set<String>::const_iterator it1 = new_sequences.begin(); it1 != new_sequences.end(); ++it1)
			{
				cerr << tabs_ << *it1 << endl;
			}
			cerr << endl;
#endif
		}

		for (set<String>::const_iterator sit = new_sequences.begin(); sit != new_sequences.end(); ++sit)
		{
			sequences.insert(*sit);
		}
  }
#ifdef DAC_DEBUG
	cerr << tabs_ << "Found sequences for " << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << endl;
  for (set<String>::const_iterator sit = sequences.begin(); sit != sequences.end(); ++sit)
	{
		cerr << tabs_ << *sit << endl;
	}
#endif
	
	// reduce the permuts once again to reduce complexity
	if (offset_prefix > precursor_mass_tolerance_ || offset_suffix > precursor_mass_tolerance_)
	{
		reducePermuts_(sequences, CID_spec, offset_prefix, offset_suffix);
	}

#ifdef DAC_DEBUG
	cerr << "Writing to cache " << left << " " << right << endl;
#endif
	
	subspec_to_sequences_[left][right] = sequences;

#ifdef DAC_DEBUG
	depth_--;
	cerr << tabs_ << "DAC: " << CID_spec[left].getPosition()[0] << " " << CID_spec[right].getPosition()[0] << " " << sequences.size() << endl;
#endif
	return;

}

	void CompNovoIdentificationCID::updateMembers_()
	{
		CompNovoIdentificationBase::updateMembers_();
		precursor_mass_tolerance_ = param_.getValue("precursor_mass_tolerance");
		return;
	}

}
