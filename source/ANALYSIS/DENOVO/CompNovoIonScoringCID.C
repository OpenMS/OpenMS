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

#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringCID.h>
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecompositionAlgorithm.h>
#include <OpenMS/CONCEPT/Constants.h>

//#define ION_SCORING_DEBUG
//#define SCORE_WITNESSSET_DEBUG

using namespace std;

namespace OpenMS
{
	CompNovoIonScoringCID::CompNovoIonScoringCID()
		:	CompNovoIonScoringBase()
	{
		defaults_.setValue("precursor_mass_tolerance", 1.5, "precursor mass tolerance");

		defaultsToParam_();
		updateMembers_();
	}

	CompNovoIonScoringCID::CompNovoIonScoringCID(const CompNovoIonScoringCID& rhs)
		: CompNovoIonScoringBase(rhs)
	{
	}

	CompNovoIonScoringCID& CompNovoIonScoringCID::operator = (const CompNovoIonScoringCID& rhs)
	{
		if (this != &rhs)
		{
			CompNovoIonScoringBase::operator = (rhs);
		}
		return *this;
	}

	CompNovoIonScoringCID::~CompNovoIonScoringCID()
	{
	}

	void CompNovoIonScoringCID::scoreSpectrum(Map<DoubleReal, IonScore>& ion_scores, PeakSpectrum& CID_spec, DoubleReal precursor_weight, Size charge)
	{
    for (PeakSpectrum::ConstIterator it = CID_spec.begin(); it != CID_spec.end(); ++it)
    {
      DoubleReal it_pos(it->getPosition()[0]);
			IonScore ion_score;
			ion_scores[it_pos] = ion_score;
    }

		// adds single charged variants of putative single charged ions
		//addSingleChargedIons_(ion_scores, CID_spec);

		for (PeakSpectrum::ConstIterator it = CID_spec.begin(); it != CID_spec.end(); ++it)
		{
			DoubleReal it_pos(it->getPosition()[0]);
			IonScore ion_score;
			ion_scores[it_pos] = ion_score;
		}		

    for (PeakSpectrum::ConstIterator it = CID_spec.begin(); it != CID_spec.end(); ++it)
    {
      ion_scores[it->getPosition()[0]].s_isotope_pattern_1 = scoreIsotopes_(CID_spec, it, ion_scores, 1);
      if (it->getPosition()[0] < precursor_weight / 2.0)
      {
        ion_scores[it->getPosition()[0]].s_isotope_pattern_2 =  scoreIsotopes_(CID_spec, it, ion_scores, 2);
      }
      else
      {
        ion_scores[it->getPosition()[0]].s_isotope_pattern_2 = -1;
      }
    }

    // combine the features and give y-ion scores
    scoreWitnessSet_(charge, precursor_weight, ion_scores, CID_spec);

    for (Map<DoubleReal, IonScore>::iterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
    {
      it->second.score = it->second.s_witness;
    }

		
		MassDecompositionAlgorithm decomp_algo;
		Param decomp_param(decomp_algo.getParameters());
		decomp_param.setValue("tolerance", fragment_mass_tolerance_);
		decomp_algo.setParameters(decomp_param);

		DoubleReal y_offset = EmpiricalFormula("H2O").getMonoWeight() + Constants::PROTON_MASS_U;
		
		// check whether a PRMNode_ can be decomposed into amino acids
    // rescore the peaks that cannot be possible y-ion candidates
		DoubleReal max_decomp_weight(param_.getValue("max_decomp_weight"));
    for (Map<DoubleReal, IonScore>::iterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
    {
      if (it->first > y_offset && (it->first - y_offset) < max_decomp_weight)
      {
        vector<MassDecomposition> decomps;
        decomp_algo.getDecompositions(decomps, it->first - y_offset);
#ifdef ION_SCORING_DEBUG
        cerr << "Decomps: " << it->first <<  " " << it->first - y_offset << " " << decomps.size() << " " << it->second.score << endl;
#endif
        if (decomps.empty())
        {
          it->second.score = 0;
        }
      }
		}

		decomp_param.setValue("tolerance", (DoubleReal)param_.getValue("precursor_mass_tolerance"));
		decomp_algo.setParameters(decomp_param);
		// now the upper part with differen tolerance
		for (Map<DoubleReal, IonScore>::iterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
		{
      if (it->first < precursor_weight && precursor_weight - it->first < max_decomp_weight)
      {
        vector<MassDecomposition> decomps;
        decomp_algo.getDecompositions(decomps, precursor_weight - it->first);
#ifdef ION_SCORING_DEBUG
        cerr << "Decomps: " << it->first << " " << precursor_weight - it->first << " " << decomps.size() << " " << it->second.score << endl;
#endif
        if (decomps.empty())
        {
          it->second.score = 0;
        }
      }
    }

    ion_scores[CID_spec.begin()->getPosition()[0]].score = 1;
    ion_scores[(CID_spec.end() - 1)->getPosition()[0]].score = 1;
	}

void CompNovoIonScoringCID::scoreWitnessSet_(Size charge, DoubleReal precursor_weight, Map<DoubleReal, IonScore>& ion_scores, const PeakSpectrum& CID_spec)
{
	DoubleReal precursor_mass_tolerance = (DoubleReal)param_.getValue("precursor_mass_tolerance");
	vector<DoubleReal> diffs;
	//diffs.push_back(28.0);
	diffs.push_back(EmpiricalFormula("NH3").getMonoWeight());
	diffs.push_back(EmpiricalFormula("H2O").getMonoWeight());

	// witnesses of CID spec (diffs)
	for (PeakSpectrum::ConstIterator it1 = CID_spec.begin(); it1 != CID_spec.end(); ++it1)
	{
		//Size num_wit(0);
		DoubleReal wit_score(0.0);
		DoubleReal pos1(it1->getPosition()[0]);
		wit_score += it1->getIntensity();
		for (PeakSpectrum::ConstIterator it2 = CID_spec.begin(); it2 != CID_spec.end(); ++it2)
		{
			DoubleReal pos2(it2->getPosition()[0]);

			// direct ++
			if (charge > 1)
			{
				if (fabs(pos2 * 2 - Constants::PROTON_MASS_U - pos1) < fragment_mass_tolerance_)
				{
					DoubleReal factor((fragment_mass_tolerance_ - fabs(pos2 * 2 - Constants::PROTON_MASS_U - pos1)) / fragment_mass_tolerance_);
					// pos1 is ion, pos2 is ++ion
#ifdef SCORE_WITNESSSET_DEBUG
					cerr << "scoreWitnessSet: ++ion " << pos1 << " " << pos2 << " (factor=" << factor << ") " << wit_score << " -> ";
#endif
        	if (ion_scores[it2->getPosition()[0]].s_isotope_pattern_2 < 0.2)
        	{
          	wit_score += it2->getIntensity() * 0.2 * factor;
        	}
        	else
        	{
          	wit_score += it2->getIntensity() * ion_scores[it2->getPosition()[0]].s_isotope_pattern_2 * factor;
        	}
#ifdef SCORE_WITNESSSET_DEBUG
					cerr << wit_score << endl;
#endif
				}
			}

			// diffs?
			for (vector<DoubleReal>::const_iterator it = diffs.begin(); it != diffs.end(); ++it)
			{
				// pos1 is ion, pos2 loss peak
				if (fabs(pos1 - pos2 - *it) < precursor_mass_tolerance)
				{
					DoubleReal factor((fragment_mass_tolerance_ - fabs(pos1 - pos2 - *it)) / fragment_mass_tolerance_);
#ifdef SCORE_WITNESSSET_DEBUG
					cerr << "scoreWitnessSet: diff " << pos1 << " (" << pos2 << ") " << *it << " (factor=" << factor << ") " << wit_score << " -> ";
#endif
					wit_score += it2->getIntensity() / 5.0 * factor;
#ifdef SCORE_WITNESSSET_DEBUG
					cerr << wit_score << endl;
#endif
				}
			}

			// is there a b-ion?; pos1 is ion, pos2 complementary ion
			if (fabs(pos1 + pos2 - 1 * Constants::PROTON_MASS_U - precursor_weight) < fragment_mass_tolerance_)
			{
				DoubleReal factor((fragment_mass_tolerance_ - fabs(pos1 + pos2 - Constants::PROTON_MASS_U - precursor_weight)) / fragment_mass_tolerance_);
				factor *= 0.2;
#ifdef SCORE_WITNESSSET_DEBUG
				cerr << "scoreWitnessSet: complementary " << pos1 << " (" << pos2 << ") (factor=" << factor << ") " << wit_score << " -> ";
#endif
				// found complementary ion
				if (ion_scores[it2->getPosition()[0]].s_isotope_pattern_1 < 0.5 || ion_scores[it2->getPosition()[0]].is_isotope_1_mono != 1)
				{
					wit_score += it2->getIntensity() * 0.5 * factor;
				}
				else
				{
					wit_score += it2->getIntensity() * ion_scores[it2->getPosition()[0]].s_isotope_pattern_1 * factor;
				}
#ifdef SCORE_WITNESSSET_DEBUG
				cerr << wit_score << endl;
#endif

				if (ion_scores[it2->getPosition()[0]].s_bion != 0)
				{
#ifdef SCORE_WITNESSSET_DEBUG
					cerr << "scoreWitnessSet: complementary is b-ion " << pos1 << "(" << pos2 << ")" << wit_score << " -> ";
#endif
					wit_score += ion_scores[it2->getPosition()[0]].s_bion * factor;
#ifdef SCORE_WITNESSSET_DEBUG
					cerr << wit_score << endl;
#endif
				}

			}
		}

		// isotope pattern ok?
		if (ion_scores[it1->getPosition()[0]].s_isotope_pattern_1 > 0 && ion_scores[it1->getPosition()[0]].is_isotope_1_mono == 1)
		{
#ifdef SCORE_WITNESSSET_DEBUG
			cerr << "scoreWitnessSet: isotope pattern: " << pos1 << " " << wit_score << " -> ";
#endif
			wit_score += ion_scores[it1->getPosition()[0]].s_isotope_pattern_1 * wit_score;
#ifdef SCORE_WITNESSSET_DEBUG
			cerr << wit_score << endl;
#endif
		}
		
		if (ion_scores[it1->getPosition()[0]].s_yion > 0)
		{
#ifdef SCORE_WITNESSSET_DEBUG
			cerr << "scoreWitnessSet: is y-ion: " << pos1 << " " << wit_score << " -> ";
#endif
			wit_score += ion_scores[it1->getPosition()[0]].s_yion;
#ifdef SCORE_WITNESSSET_DEBUG
			cerr << wit_score << endl;
#endif
		}

		if (ion_scores[it1->getPosition()[0]].s_bion > 0)
		{
#ifdef SCORE_WITNESSSET_DEBUG
			cerr << "scoreWitnessSet: is b-ion: " << pos1 << " " << wit_score << " -> ";
#endif
			if (ion_scores[it1->getPosition()[0]].s_bion < wit_score)
			{
				wit_score -= ion_scores[it1->getPosition()[0]].s_bion;
			}
			else
			{
				wit_score = 0;
			}
		}
		
		ion_scores[it1->getPosition()[0]].s_witness = wit_score;
	}
	return;
}

}


