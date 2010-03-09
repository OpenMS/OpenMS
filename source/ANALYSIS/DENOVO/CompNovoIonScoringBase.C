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

#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoringBase.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <numeric>

//#define ION_SCORING_DEBUG

using namespace std;

namespace OpenMS
{
	CompNovoIonScoringBase::IonScore::IonScore()
		: score(0),
      s_bion(0),
      s_yion(0),
      s_witness(0),
      position(0),
      s_isotope_pattern_1(0),
      is_isotope_1_mono(0),
      s_isotope_pattern_2(0)
	{
	}

	CompNovoIonScoringBase::IonScore::IonScore(const IonScore& rhs)
    : score(rhs.score),
      s_bion(rhs.s_bion),
      s_yion(rhs.s_yion),
      s_witness(rhs.s_witness),
      position(rhs.position),
      s_isotope_pattern_1(rhs.s_isotope_pattern_1),
      is_isotope_1_mono(rhs.is_isotope_1_mono),
      s_isotope_pattern_2(rhs.s_isotope_pattern_2)
	{
	}

	CompNovoIonScoringBase::IonScore::~IonScore()
	{
	}

	CompNovoIonScoringBase::IonScore& CompNovoIonScoringBase::IonScore::operator = (const IonScore& rhs)
	{
		if (this != &rhs)
		{
			score = rhs.score;
			s_bion = rhs.s_bion;
			s_yion = rhs.s_yion;
			s_witness = rhs.s_witness;
			position = rhs.position;
			s_isotope_pattern_1 = rhs.s_isotope_pattern_1;
			is_isotope_1_mono = rhs.is_isotope_1_mono;
			s_isotope_pattern_2 = rhs.s_isotope_pattern_2;
		}
		return *this;
	}
				
	CompNovoIonScoringBase::CompNovoIonScoringBase()
		:	DefaultParamHandler("CompNovoIonScoringBase"),
     	fragment_mass_tolerance_(0)
	{
		defaults_.setValue("fragment_mass_tolerance", 0.4, "fragment mass tolerance");
		defaults_.setValue("decomp_weights_precision", 0.01, "precision used to calculate the decompositions, this only affects cache usage!", StringList::create("advanced"));
		defaults_.setValue("double_charged_iso_threshold", 0.9, "minimal isotope intensity correlation of doubly charged ions to be used to score the single scored ions", StringList::create("advanced"));
		defaults_.setValue("double_charged_iso_threshold_single", 0.99, "Isotope scoring threshold used for doubly charged ions to infer singly charged variants", StringList::create("advanced"));
		defaults_.setValue("max_isotope_to_score", 3, "max isotope peak to be considered in the scoring", StringList::create("advanced"));
		defaults_.setValue("max_decomp_weight", 600, "maximal m/z difference used to calculate the decompositions", StringList::create("advanced"));
		defaults_.setValue("max_isotope", 3, "max isotope used in the theoretical spectra to score", StringList::create("advanced"));
		defaults_.setValue("max_mz", 2000.0, "maximal m/z value used to calculate isotope distributions", StringList::create("advanced"));

		defaultsToParam_();
	}

	CompNovoIonScoringBase::CompNovoIonScoringBase(const CompNovoIonScoringBase& rhs)
		: DefaultParamHandler(rhs)
	{
		updateMembers_();
	}

	CompNovoIonScoringBase& CompNovoIonScoringBase::operator = (const CompNovoIonScoringBase& rhs)
	{
		if (this != &rhs)
		{
			DefaultParamHandler::operator=(rhs);
			updateMembers_();
			// TODO
		}
		return *this;
	}

	CompNovoIonScoringBase::~CompNovoIonScoringBase()
	{
	}

	void CompNovoIonScoringBase::addSingleChargedIons_(Map<DoubleReal, IonScore>& ion_scores, PeakSpectrum& CID_spec)
	{
		DoubleReal double_charged_iso_threshold_single((DoubleReal)param_.getValue("double_charged_iso_threshold_single"));
		PeakSpectrum CID_spec_new = CID_spec;
		for (PeakSpectrum::ConstIterator it = CID_spec.begin(); it != CID_spec.end(); ++it)
		{
			if (it->getPosition()[0] < CID_spec.getPrecursors().begin()->getMZ() / 2.0)	
			{
				DoubleReal score = scoreIsotopes_(CID_spec, it, ion_scores, 2);
				if (score > double_charged_iso_threshold_single)
				{
					// infer this peak as single charged variant
					DoubleReal mz_comp = it->getPosition()[0] * 2.0 - Constants::PROTON_MASS_U;
					bool found(false);
					for (PeakSpectrum::ConstIterator it1 = CID_spec.begin(); it1 != CID_spec.end(); ++it1)
					{
						if (fabs(mz_comp - it1->getPosition()[0]) < fragment_mass_tolerance_)
						{
							found = true;
							break;
						}
					}

					if (!found)
					{	
						Peak1D p;
						p.setIntensity(it->getIntensity());
						p.setPosition(mz_comp);
						CID_spec_new.push_back(p);
					}
				}
			}
			else
			{
				break;
			}
		}

		CID_spec = CID_spec_new;
	}
	
	CompNovoIonScoringBase::IsotopeType CompNovoIonScoringBase::classifyIsotopes_(const PeakSpectrum& spec, PeakSpectrum::ConstIterator it)
	{
  	DoubleReal it_pos(it->getPosition()[0]);

  	// is there a peak left of it with diff 1Da?
  	for (PeakSpectrum::ConstIterator it1 = it; it1 != spec.end(); --it1)
  	{
    	DoubleReal it1_pos(it1->getPosition()[0]);

    	if (it1 == spec.begin() || fabs(it_pos - it1_pos) > 1.5)
    	{
      	break;
    	}

    	if (fabs(fabs(it_pos - it1_pos) - 1.0) < fragment_mass_tolerance_)
    	{
    	  return CHILD;
    	}
  	}

  	// is there a peak right of it with diff 1Da?
  	for (PeakSpectrum::ConstIterator it1 = it; it1 != spec.end(); ++it1)
  	{
    	DoubleReal it1_pos(it1->getPosition()[0]);
    	if (fabs(fabs(it_pos - it1_pos) - 1.0) < fragment_mass_tolerance_)
    	{		
      	return PARENT;
    	}

    	if (fabs(it_pos - it1_pos) > 1.5)
    	{
      	break;
    	}
  	}

  	return LONE;
	}
	
	DoubleReal CompNovoIonScoringBase::scoreIsotopes_(const PeakSpectrum& CID_spec, PeakSpectrum::ConstIterator it, Map<DoubleReal, IonScore>& ion_scores, Size charge)
	{
  	DoubleReal it_pos(it->getMZ());  // ~ weight of the fragment
		UInt max_isotope_to_score(param_.getValue("max_isotope_to_score"));
		DoubleReal double_charged_iso_threshold(param_.getValue("double_charged_iso_threshold"));
  	DoubleReal actual_pos = it_pos;

  	vector<DoubleReal> iso_pattern;
  	vector<PeakSpectrum::ConstIterator> iso_pattern_its;
  	iso_pattern.push_back(it->getIntensity());
  	iso_pattern_its.push_back(it);
  	// get all peaks that have the right distance right of the given peak
  	for (PeakSpectrum::ConstIterator it1 = it; it1 != CID_spec.end(); ++it1)
  	{
    	DoubleReal it1_pos(it1->getPosition()[0]);
    	if (fabs(fabs(actual_pos - it1_pos) - Constants::NEUTRON_MASS_U / (DoubleReal)charge) < fragment_mass_tolerance_)
    	{
      	iso_pattern.push_back(it1->getIntensity());
      	actual_pos = it1_pos;
      	iso_pattern_its.push_back(it1);
    	}

    	if (iso_pattern.size() == max_isotope_to_score)
    	{
      	break;
    	}
  	}

  	if (iso_pattern.size() == 1)
  	{
    	return -1;
  	}

  	// normalize the intensity to a sum of one
  	DoubleReal sum(0);
  	for (vector<DoubleReal>::const_iterator it1 = iso_pattern.begin(); it1 != iso_pattern.end(); ++it1)
  	{
    	sum += *it1;
  	}

  	for (vector<DoubleReal>::iterator it1 = iso_pattern.begin(); it1 != iso_pattern.end(); ++it1)
  	{
    	*it1 /= sum;
  	}

  	// get the theoretical isotope distribution
  	IsotopeDistribution iso_dist(iso_pattern.size());
  	iso_dist.estimateFromPeptideWeight((it_pos - charge * Constants::PROTON_MASS_U) * charge + Constants::PROTON_MASS_U);

  	// compare the distribution sizes
  	if (iso_dist.size() != iso_pattern.size())
  	{
    	cerr << "scoreIsotopes: error istope distributions have differing sizes" << endl;
    	return -1;
  	}

	  // calculate simple correlation score
  	DoubleReal score(0.0);

  	DoubleReal numerator(0), auto1(0), auto2(0);
  	for (Size i = 0; i != iso_dist.size(); ++i)
  	{
    	numerator += iso_dist.getContainer()[i].second * iso_pattern[i];
    	auto1 += iso_dist.getContainer()[i].second * iso_dist.getContainer()[i].second;
    	auto2 += iso_pattern[i] * iso_pattern[i];
  	}

  	score = numerator * numerator / auto1 / auto2;

  	// if score is great enough, we accept it
  	if (score > double_charged_iso_threshold)
  	{
    	if (ion_scores[it_pos].is_isotope_1_mono == 0)
    	{
      	ion_scores[it_pos].is_isotope_1_mono = 1;
    	}

    	for (Size i = 1; i < iso_pattern_its.size(); ++i)
    	{
      	ion_scores[iso_pattern_its[i]->getPosition()[0]].is_isotope_1_mono = -1;
#ifdef ION_SCORING_DEBUG
      	cerr << "scoreIsotopes: disabling " << iso_pattern_its[i]->getPosition()[0] << endl;
#endif
    	}
  	}
#ifdef ION_SCORING_DEBUG
  	cerr << "IsotopeScore: " << it_pos << " " << score << " " << iso_dist.size() << " " << ion_scores[it->getPosition()[0]].is_isotope_1_mono << " z=" << charge << endl;
#endif
  	return score;
	}


DoubleReal CompNovoIonScoringBase::scoreIsotopes(const PeakSpectrum& spec, PeakSpectrum::ConstIterator it, Size charge)
{
#ifdef ION_SCORING_DEBUG
	cerr << "scoreIsotopes: " << spec.size() << " " << it->getPosition()[0] << " " << it->getIntensity() << " " << charge << endl;
#endif
  DoubleReal it_pos(it->getMZ()); // ~ weight of the fragment
  DoubleReal actual_pos = it_pos;
	UInt max_isotope_to_score = (UInt)param_.getValue("max_isotope_to_score");

  vector<DoubleReal> iso_pattern;
  iso_pattern.push_back(it->getIntensity());

	// get all peaks that have the right distance right of the given peak
	//cerr << "Scoring peaks...";
  for (PeakSpectrum::ConstIterator it1 = it; it1 != spec.end(); ++it1)
  {
    DoubleReal it1_pos(it1->getMZ());

		//cerr << "PRE: " << actual_pos << " " << it1_pos << " " << Constants::NEUTRON_MASS_U << " " << charge << " " << fragment_mass_tolerance_ << endl;
		
    if (fabs(fabs(actual_pos - it1_pos) - Constants::NEUTRON_MASS_U / (DoubleReal)charge) < fragment_mass_tolerance_ / (DoubleReal)charge)
    {
#ifdef ION_SCORING_DEBUG
			cerr << actual_pos << " " << it1_pos << " " << charge << " " << fragment_mass_tolerance_ << endl;
#endif
      iso_pattern.push_back(it1->getIntensity());
      actual_pos = it1_pos;
    }

    if (iso_pattern.size() == max_isotope_to_score)
    {
      break;
    }
  }
	//cerr << "ended" << endl;

  if (iso_pattern.size() == 1)
  {
    return 0;
  }

  // normalize the intensity to a sum of one
	/*
  DoubleReal sum(0);
  for (vector<DoubleReal>::const_iterator it1 = iso_pattern.begin(); it1 != iso_pattern.end(); ++it1)
  {
    sum += *it1;
  }

  for (vector<DoubleReal>::iterator it1 = iso_pattern.begin(); it1 != iso_pattern.end(); ++it1)
  {
    *it1 /= sum;
  }*/


  // get the theoretical isotope distribution
  IsotopeDistribution iso_dist(iso_pattern.size());
  iso_dist.estimateFromPeptideWeight(it_pos * (DoubleReal)charge - (DoubleReal)(charge - 1) * Constants::PROTON_MASS_U);

  // compare the distribution sizes
  if (iso_dist.size() != iso_pattern.size())
  {
    cerr << "scoreIsotopes: error istope distributions have differing sizes" << endl;
    return -1;
  }

  // calculate simple correlation score
  DoubleReal score(0.0);

  DoubleReal numerator(0), auto1(0), auto2(0);
  for (Size i = 0; i != iso_dist.size(); ++i)
  {
    numerator += iso_dist.getContainer()[i].second * iso_pattern[i];
    auto1 += iso_dist.getContainer()[i].second * iso_dist.getContainer()[i].second;
    auto2 += iso_pattern[i] * iso_pattern[i];
  }

	//score *= accumulate(iso_pattern.begin(), iso_pattern.end(), 0.0);

  score = numerator * numerator / auto1 / auto2;
	score *= accumulate(iso_pattern.begin(), iso_pattern.end(), 0.0);
#ifdef ION_SCORING_DEBUG
  cerr << "IsotopeScore: " << it_pos << " " << score << " " << iso_dist.size() << " z=" << charge << endl;
#endif
  return score;
}

	void CompNovoIonScoringBase::initIsotopeDistributions_()
	{
		DoubleReal max_mz(param_.getValue("max_mz"));
		UInt max_isotope(param_.getValue("max_isotope"));
  	IsotopeDistribution iso_dist(max_isotope);
  	for (Size i = 1; i <= max_mz; ++i)
  	{
    iso_dist.estimateFromPeptideWeight((DoubleReal)i);
    iso_dist.renormalize();
    vector<DoubleReal> iso(max_isotope, 0.0);

    for (Size j = 0; j != iso_dist.size(); ++j)
    {
      iso[j] = iso_dist.getContainer()[j].second;
    }
    	isotope_distributions_[i] = iso;
  	}
	}

	void 	CompNovoIonScoringBase::updateMembers_()
	{
		fragment_mass_tolerance_ = (DoubleReal)param_.getValue("fragment_mass_tolerance");

		initIsotopeDistributions_();

		return;
	}
}


