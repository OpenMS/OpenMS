// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoring.h>
#include <OpenMS/ANALYSIS/DENOVO/MassDecompositionAlgorithm.h>
#include <OpenMS/CONCEPT/Constants.h>

//#define ION_SCORING_DEBUG
//#define SCORE_ETDFEATURES_DEBUG
//#define SCORE_WITNESSSET_DEBUG

using namespace std;

namespace OpenMS
{
	CompNovoIonScoring::CompNovoIonScoring()
		:	CompNovoIonScoringBase()
	{
	}

	CompNovoIonScoring::CompNovoIonScoring(const CompNovoIonScoring& rhs)
		: CompNovoIonScoringBase(rhs)
	{
	}

	CompNovoIonScoring& CompNovoIonScoring::operator = (const CompNovoIonScoring& rhs)
	{
		if (this != &rhs)
		{
			CompNovoIonScoringBase::operator = (rhs);
		}
		return *this;
	}

	CompNovoIonScoring::~CompNovoIonScoring()
	{
	}

	void CompNovoIonScoring::scoreSpectra(Map<DoubleReal, IonScore>& ion_scores, PeakSpectrum& CID_spec, PeakSpectrum& ETD_spec, DoubleReal precursor_weight, Size charge)
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

    // find possible supporting ions from ETD spec to CID spec
    scoreETDFeatures_(charge, precursor_weight, ion_scores, CID_spec, ETD_spec);

    // combine the features and give b-ion scores
    scoreWitnessSet_(charge, precursor_weight, ion_scores, CID_spec);

    for (Map<DoubleReal, IonScore>::iterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
    {
      it->second.score = it->second.s_witness;
    }

		
		MassDecompositionAlgorithm decomp_algo;

		
		// check whether a PRMNode_ can be decomposed into amino acids
    // rescore the peaks that cannot be possible y-ion candidates
		DoubleReal max_decomp_weight((DoubleReal)param_.getValue("max_decomp_weight"));
    for (Map<DoubleReal, IonScore>::iterator it = ion_scores.begin(); it != ion_scores.end(); ++it)
    {
      if (it->first > 19.0 && (it->first - 19.0) < max_decomp_weight)
      {
        vector<MassDecomposition> decomps;
        decomp_algo.getDecompositions(decomps, it->first - 19.0);
#ifdef ION_SCORING_DEBUG
        cerr << "Decomps: " << it->first <<  " " << it->first -19.0 << " " << decomps.size() << " " << it->second.score << endl;
#endif
        if (decomps.size() == 0)
        {
          it->second.score = 0;
        }
      }

      if (it->first < precursor_weight && precursor_weight - it->first < max_decomp_weight)
      {
        vector<MassDecomposition> decomps;
        decomp_algo.getDecompositions(decomps, precursor_weight - it->first);
#ifdef ION_SCORING_DEBUG
        cerr << "Decomps: " << it->first << " " << precursor_weight - it->first << " " << decomps.size() << " " << it->second.score << endl;
#endif
        if (decomps.size() == 0)
        {
          it->second.score = 0;
        }
      }
    }

    ion_scores[CID_spec.begin()->getPosition()[0]].score = 1;
    ion_scores[(CID_spec.end() - 1)->getPosition()[0]].score = 1;
	}

void CompNovoIonScoring::scoreWitnessSet_(Size charge, DoubleReal precursor_weight, Map<DoubleReal, IonScore>& ion_scores, const PeakSpectrum& CID_spec)
{
	vector<DoubleReal> diffs;
	//diffs.push_back(28.0);
	diffs.push_back(17.0);
	diffs.push_back(18.0);

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
          	wit_score += it2->getIntensity() * /* 0.2 */ factor;
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
				if (fabs(pos1 - pos2 - *it) < fragment_mass_tolerance_)
				{
					DoubleReal factor((fragment_mass_tolerance_ - fabs(pos1 - pos2 - *it)) / fragment_mass_tolerance_);
#ifdef SCORE_WITNESSSET_DEBUG
					cerr << "scoreWitnessSet: diff " << pos1 << " (" << pos2 << ") " << *it << " (factor=" << factor << ") " << wit_score << " -> ";
#endif
					wit_score += it2->getIntensity()/* / 5.0*/ * factor;
#ifdef SCORE_WITNESSSET_DEBUG
					cerr << wit_score << endl;
#endif
				}
			}

			// is there a b-ion?; pos1 is ion, pos2 complementary ion
			if (fabs(pos1 + pos2 - 1 * Constants::PROTON_MASS_U - precursor_weight) < fragment_mass_tolerance_)
			{
				DoubleReal factor((fragment_mass_tolerance_ - fabs(pos1 + pos2 - Constants::PROTON_MASS_U - precursor_weight)) / fragment_mass_tolerance_);
				/*factor *= 0.2;*/
#ifdef SCORE_WITNESSSET_DEBUG
				cerr << "scoreWitnessSet: complementary " << pos1 << " (" << pos2 << ") (factor=" << factor << ") " << wit_score << " -> ";
#endif
				// found complementary ion
				if (ion_scores[it2->getPosition()[0]].s_isotope_pattern_1 < 0.5 || ion_scores[it2->getPosition()[0]].is_isotope_1_mono != 1)
				{
					wit_score += it2->getIntensity() /* * 0.5*/ * factor;
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


void CompNovoIonScoring::scoreETDFeatures_(Size /*charge*/, DoubleReal precursor_weight, Map<DoubleReal, IonScore>& ion_scores, const PeakSpectrum& CID_spec, const PeakSpectrum& ETD_spec)
{
	//DoubleReal fragment_mass_tolerance((DoubleReal)param_.getValue("fragment_mass_tolerance"));
	Size max_isotope_to_score(param_.getValue("max_isotope_to_score"));

	for (PeakSpectrum::ConstIterator it1 = CID_spec.begin(); it1 != CID_spec.end(); ++it1)
	{
		DoubleReal pos1(it1->getPosition()[0]);
		DoubleReal b_sum(0.0), y_sum(0.0);

		// score a-ions
		for (PeakSpectrum::ConstIterator it2 = CID_spec.begin(); it2 != CID_spec.end(); ++it2)
		{
			DoubleReal pos2(it2->getPosition()[0]);
			if (fabs(pos1 - pos2 - 28.0) < fragment_mass_tolerance_)
			{
				DoubleReal factor((fragment_mass_tolerance_ - fabs(pos1 - pos2 - 28.0)) / fragment_mass_tolerance_);
#ifdef SCORE_ETDFEATURES_DEBUG
				cerr << "scoreETDFeatures: found a-ion " << pos1 << " (" << pos2 << ") (factor=" << factor << ") " << b_sum << " -> ";
#endif
				b_sum += it2->getIntensity() * factor;
#ifdef SCORE_ETDFEATURES_DEBUG
				cerr << endl;
#endif
			}
		}
		
		for (PeakSpectrum::ConstIterator it2 = ETD_spec.begin(); it2 != ETD_spec.end(); ++it2)
		{
			DoubleReal pos2(it2->getPosition()[0]);

			// check if pos2 is precursor doubly charged, which has not fragmented
			DoubleReal pre_diff_lower = (precursor_weight + Constants::PROTON_MASS_U) / 2.0 - fragment_mass_tolerance_;
			DoubleReal pre_diff_upper = (precursor_weight + 4.0 * Constants::PROTON_MASS_U) / 2.0 + fragment_mass_tolerance_;
			if (pos2 > pre_diff_lower && pos2 < pre_diff_upper)
			{
#ifdef SCORE_ETDFEATURES_DEBUG
				cerr << "scoreETDFeatures: pre-range: " << pos2 << " is in precursor peak range: " << pre_diff_lower << " <-> " << pre_diff_upper << endl;
#endif
				continue;
			}
			
			//DoubleReal diff(pos2 - pos1);

			// pos1 is CID ion; pos2 is ETD ion
			// pos1 b-ion, pos2 c-ion
			if (fabs(pos1 + 17.0 - pos2) < fragment_mass_tolerance_)
			{
				// now test if the ETD peak has "isotope" pattern
				DoubleReal factor((fragment_mass_tolerance_ - fabs(pos1 + 17.0 - pos2)) / fragment_mass_tolerance_);
#ifdef SCORE_ETDFEATURES_DEBUG
				cerr << "scoreETDFeatures: is b-ion: " << pos1 << " (" << pos2 << ") (factor=" << factor << ") " << b_sum << " -> ";
#endif
				vector<DoubleReal> iso_pattern;
				iso_pattern.push_back(it1->getIntensity());
				DoubleReal actual_pos = it1->getPosition()[0];
				for (PeakSpectrum::ConstIterator it3 = it2; it3 != ETD_spec.end(); ++it3)
				{
					DoubleReal it3_pos(it3->getPosition()[0]);
					if (fabs(fabs(actual_pos - it3_pos) - Constants::NEUTRON_MASS_U) < fragment_mass_tolerance_)
					{
						iso_pattern.push_back(it3->getIntensity());
						actual_pos = it3_pos;
					}
					if (iso_pattern.size() == max_isotope_to_score)
					{
						break;
					}
				}
				
				if (ion_scores[it1->getPosition()[0]].is_isotope_1_mono != -1)
				{
					b_sum += it2->getIntensity() * iso_pattern.size() * factor;
				}
#ifdef SCORE_ETDFEATURES_DEBUG					
				cerr << b_sum << endl;
#endif
			}

			
			
			// pos1 z-ion, pos2 y-ion
			if (fabs(pos2 + 16.0 - pos1) < fragment_mass_tolerance_)
			{
				DoubleReal factor((fragment_mass_tolerance_ - fabs(pos2 + 16.0 - pos1)) / fragment_mass_tolerance_);
				// now test if the ETD peak has "isotope" pattern
#ifdef SCORE_ETDFEATURES_DEBUG
				cerr << "scoreETDFeatures: is y-ion: " << pos1 << " (" << pos2 << ") (factor=" << factor << ") " << y_sum << " -> ";
#endif
        vector<DoubleReal> iso_pattern;
        iso_pattern.push_back(it1->getIntensity());
        DoubleReal actual_pos = it1->getPosition()[0];
        for (PeakSpectrum::ConstIterator it3 = it2; it3 != ETD_spec.end(); ++it3)
        {
          DoubleReal it3_pos(it3->getPosition()[0]);
          if (fabs(fabs(actual_pos - it3_pos) - Constants::NEUTRON_MASS_U) < fragment_mass_tolerance_)
          {
            iso_pattern.push_back(it3->getIntensity());
            actual_pos = it3_pos;
          }
          if (iso_pattern.size() == max_isotope_to_score)
          {
            break;
          }
        }
#ifdef SCORE_ETDFEATURES_DEBUG
				cerr << ion_scores[it1->getPosition()[0]].is_isotope_1_mono << " ";
#endif
				if (ion_scores[it1->getPosition()[0]].is_isotope_1_mono != -1)
				{
					y_sum += it2->getIntensity() * iso_pattern.size() * factor;
				}
#ifdef SCORE_ETDFEATURES_DEBUG
				cerr << y_sum << endl;
#endif
			}
		}
		ion_scores[it1->getPosition()[0]].s_bion = b_sum;
		ion_scores[it1->getPosition()[0]].s_yion = y_sum;
	}
	return;
}

}


