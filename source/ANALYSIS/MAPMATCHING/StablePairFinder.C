// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/StablePairFinder.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/SYSTEM/StopWatch.h>
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>

#ifdef Debug_StablePairFinder
#define V_(bla) std::cout << __FILE__ ":" << __LINE__ << ": " << bla << std::endl;
#else
#define V_(bla) ;
#endif
#define VV_(bla) V_(""#bla": " << bla);

using namespace std;

namespace OpenMS
{

  StablePairFinder::StablePairFinder() :
    Base()
  {
    //set the name for DefaultParamHandler error messages
    Base::setName(getProductName());

    defaults_.setValue("diff_exponent:RT", 1.0,
      "RT differences are raised to this power", StringList::create("advanced"));
    defaults_.setMinFloat("diff_exponent:RT",0.1);

    defaults_.setValue("diff_exponent:MZ", 2.0,
      "MZ differences are raised to this power", StringList::create("advanced"));
    defaults_.setMinFloat("diff_exponent:MZ",0.1);

    defaults_.setSectionDescription("diff_exponent",
      "Absolute position differences are raised to this power. "
        "E.g. 1 for 'linear' distance, 2 for 'quadratic' distance");

    defaults_.setValue("intensity_exponent", 0.0,
      "Intensity ratios are raised to this power.  "
        "If set to 0, intensities are not considered.", StringList::create(
        "advanced"));
    defaults_.setMinFloat("intensity_exponent",0.);

    defaults_.setValue("max_pair_distance:RT", 100.0,
      "Maximal allowed distance in RT for a pair");
    defaults_.setMinFloat("max_pair_distance:RT",0.);

    defaults_.setValue("max_pair_distance:MZ", 0.3, "Maximal allowed distance in MZ for a pair [Unit defined by 'mz_unit']");
    defaults_.setMinFloat("max_pair_distance:MZ",0.);
    defaults_.setValue("max_pair_distance:MZ_unit", "Da", "Unit of 'MZ' parameter");
    defaults_.setValidStrings("max_pair_distance:MZ_unit",StringList::create("Da,ppm"));

    defaults_.setSectionDescription("max_pair_distance",
      "Maximal allowed distance for a pair. "
        "This uses the sum of the distances in (1.) RT, (2.) MZ, and "
        "(3.) the ratio of intensities, "
        "(after modification by the corresponding "
        "diff_exponent and the intensity_exponent), and "
        "(4.) whether the charge states are equal.");

    defaults_.setValue("second_nearest_gap", 2.0,
      "The distance for the second nearest neighbors must be larger "
        "by this factor than the distance for the matching pair itself");
    defaults_.setMinFloat("second_nearest_gap",1.);

    defaults_.setValue("different_charge_penalty", 1.0,
      "If charge states are different, distance is multiplied by this factor.  "
      "If set to 1, charge has no influence.");
    defaults_.setMinFloat("different_charge_penalty",1.);

		defaults_.setValue("use_identifications", "false", "Never link features that are annotated with different peptides (only the best hit per peptide identification is taken into account).");
		defaults_.setValidStrings("use_identifications", 
															StringList::create("true,false"));

    Base::defaultsToParam_();
  }

  void
  StablePairFinder::updateMembers_()
  {
    V_("@@@ StablePairFinder::updateMembers_()");

    diff_exponent_[RT] = (DoubleReal) param_.getValue("diff_exponent:RT");
    diff_exponent_[MZ] = (DoubleReal) param_.getValue("diff_exponent:MZ");
    intensity_exponent_ = (DoubleReal) param_.getValue("intensity_exponent");
    max_pair_distance_[MZ] = (DoubleReal) param_.getValue("max_pair_distance:MZ");
    max_pair_distance_reciprocal_[MZ] = 1. / max_pair_distance_[MZ];
    max_pair_distance_mz_as_Da_ = ((String) param_.getValue("max_pair_distance:MZ_unit")) == "Da";
    max_pair_distance_[RT] = (DoubleReal) param_.getValue("max_pair_distance:RT");
    max_pair_distance_reciprocal_[RT] = 1. / max_pair_distance_[RT];
    different_charge_penalty_ = (DoubleReal) param_.getValue("different_charge_penalty");
    second_nearest_gap_ = (DoubleReal) param_.getValue("second_nearest_gap");
    return;
  }

  DoubleReal
  StablePairFinder::distance_( ConsensusFeature const & left,
                               ConsensusFeature const & right ) const
  {
    // distance from position
    DPosition<2> position_difference = left.getPosition() - right.getPosition();
    // .. in RT
    if ( position_difference[RT] < 0 )
    {
      position_difference[RT] = -position_difference[RT];
    }
    position_difference[RT] *= max_pair_distance_reciprocal_[RT];
    position_difference[RT] = pow(position_difference[RT],diff_exponent_[RT]);

    // .. in MZ
    if ( position_difference[MZ] < 0 )
    {
      position_difference[MZ] = -position_difference[MZ];
    }
    if (max_pair_distance_mz_as_Da_)
    {
      // do nothing.. we already have distance in Da
    }
    else //PPM
    {
      // distance in PPM (with 'left' as reference point)
      position_difference[MZ] /= left.getMZ();
      position_difference[MZ] *= 1e6;
    }
    position_difference[MZ] *= max_pair_distance_reciprocal_[MZ];
    position_difference[MZ] = pow(position_difference[MZ],diff_exponent_[MZ]);

    DoubleReal result = position_difference[RT] + position_difference[MZ];
    
    // distance from intensity
    if ( intensity_exponent_ != 0 )
    {
      DoubleReal right_intensity(right.getIntensity());
      DoubleReal left_intensity(left.getIntensity());
      if ( right_intensity == 0 || left_intensity == 0 )
      {
        return numeric_limits<DoubleReal>::max();
      }
      DoubleReal intensity_ratio;
      if ( left_intensity < right_intensity )
      {
        intensity_ratio = right_intensity / left_intensity;
      }
      else
      {
        intensity_ratio = left_intensity / right_intensity;
      }
      intensity_ratio = pow(intensity_ratio, intensity_exponent_);
      result *= intensity_ratio;
    }

    // distance from charge
    if ( left.getCharge() != right.getCharge() )
    {
      result *= different_charge_penalty_;
    }

    return result;
  }

  void
  StablePairFinder::run( const std::vector<ConsensusMap>& input_maps,
                         ConsensusMap &result_map )
  {
    // Empty output destination
    result_map.clear(false);

    if ( input_maps.size() != 2 )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
        "exactly two input maps required");
    }
    checkIds_(input_maps);

    std::vector<bool> is_singleton[2];
    is_singleton[0].resize(input_maps[0].size(),true);
    is_singleton[1].resize(input_maps[1].size(),true);

		bool use_IDs = String(param_.getValue("use_identifications")) == "true";

    // For each element in map 0, find its best friend in map 1
    vector<UInt> best_companion_index_0(input_maps[0].size(), UInt(-1));
    vector<DoubleReal> best_companion_distance_0(input_maps[0].size());
    vector<DoubleReal> second_best_companion_distance_0(input_maps[0].size());
    for ( UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0 )
    {
      DoubleReal best_distance = numeric_limits<DoubleReal>::max();
      DoubleReal second_best_distance = numeric_limits<DoubleReal>::max();
      for ( UInt fi1 = 0; fi1 < input_maps[1].size(); ++fi1 )
      {
        DoubleReal distance = distance_(input_maps[0][fi0], input_maps[1][fi1]);
        if ( distance < best_distance )
        {
          second_best_distance = best_distance;
          best_distance = distance;
          best_companion_index_0[fi0] = fi1;
        }
      }
      best_companion_distance_0[fi0] = best_distance;
      second_best_companion_distance_0[fi0] = second_best_distance;
    }

    // For each element in map 1, find its best friend in map 0
    vector<UInt> best_companion_index_1(input_maps[1].size(), UInt(-1));
    vector<DoubleReal> best_companion_distance_1(input_maps[1].size());
    vector<DoubleReal> second_best_companion_distance_1(input_maps[1].size());
    for ( UInt fi1 = 0; fi1 < input_maps[1].size(); ++fi1 )
    {
      DoubleReal best_distance = numeric_limits<DoubleReal>::max();
      DoubleReal second_best_distance = numeric_limits<DoubleReal>::max();
      for ( UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0 )
      {
        DoubleReal distance = distance_(input_maps[0][fi0], input_maps[1][fi1]);
        if ( distance < best_distance )
        {
          second_best_distance = best_distance;
          best_distance = distance;
          best_companion_index_1[fi1] = fi0;
        }
      }
      best_companion_distance_1[fi1] = best_distance;
      second_best_companion_distance_1[fi1] = second_best_distance;
    }

    // And if both like each other, they become a pair.
    // element_pairs_->clear();
    for ( UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0 )
    {
		// calculate separation between feature fi0 in map 0 and its best friend in map 1
		DPosition<2> position_difference = input_maps[0][fi0].getPosition() - input_maps[1][best_companion_index_0[fi0]].getPosition();
		if ( position_difference[RT] < 0 )
		{
			position_difference[RT] = -position_difference[RT];
		}
		if ( position_difference[MZ] < 0 )
		{
			position_difference[MZ] = -position_difference[MZ];
		}
		// The two features in map 0 and map 1 (fi0 and its best friend) should not be further apart then the user specified max_pair_distance.
		if ((position_difference[RT] < max_pair_distance_[RT]) && (position_difference[MZ] < max_pair_distance_[MZ]) &&
			(best_companion_distance_0[fi0] * second_nearest_gap_ <= second_best_companion_distance_0[fi0]))
		{
		// fi0 likes someone ...
        UInt best_companion_of_fi0 = best_companion_index_0[fi0];
		if ((best_companion_index_1[best_companion_of_fi0] == fi0) && 
            (best_companion_distance_1[best_companion_of_fi0] * second_nearest_gap_ <= second_best_companion_distance_1[best_companion_of_fi0]) &&
			// check if peptide IDs match:
			(!use_IDs || compatibleIDs_(input_maps[0][fi0], input_maps[1][best_companion_of_fi0])))
		{
			// ... who likes him too ...
			result_map.push_back(ConsensusFeature());
			ConsensusFeature & f = result_map.back();

			f.insert(input_maps[0][fi0]);
			f.getPeptideIdentifications().insert(f.getPeptideIdentifications().end(),input_maps[0][fi0].getPeptideIdentifications().begin(),input_maps[0][fi0].getPeptideIdentifications().end());

			f.insert(input_maps[1][best_companion_of_fi0]);
			f.getPeptideIdentifications().
            insert( f.getPeptideIdentifications().end(),input_maps[1][best_companion_of_fi0].getPeptideIdentifications().begin(),input_maps[1][best_companion_of_fi0].getPeptideIdentifications().end());

			f.computeConsensus();
			DoubleReal quality = 1. - best_companion_distance_0[fi0];
			DoubleReal quality0 = 1. - best_companion_distance_0[fi0] * second_nearest_gap_ / second_best_companion_distance_0[fi0];
			DoubleReal quality1 = 1. - best_companion_distance_1[best_companion_of_fi0] * second_nearest_gap_ / second_best_companion_distance_1[best_companion_of_fi0];
			f.setQuality(quality*quality0*quality1); // TODO other formula?

			is_singleton[0][fi0] = false;
			is_singleton[1][best_companion_of_fi0] = false;
        }
      }
    }

    // write out unmatched consensus features
    for ( UInt input = 0; input <= 1; ++ input )
    {
      for ( UInt index = 0; index < input_maps[input].size(); ++index )
      {
        if ( is_singleton[input][index] )
        {
          result_map.push_back(input_maps[input][index]);
        }
      }
    }

    // canonical ordering for checking the results, and the ids have no real meaning anyway
    result_map.sortByMZ();

		// protein IDs and unassigned peptide IDs are added to the result by the
		// FeatureGroupingAlgorithm!

    return;
  }

	
	bool StablePairFinder::compatibleIDs_(const ConsensusFeature& feat1, const ConsensusFeature& feat2) const
	{
		vector<PeptideIdentification> pep1 = feat1.getPeptideIdentifications(),pep2 = feat2.getPeptideIdentifications();
		// a feature without identifications always matches:
		if (pep1.empty() || pep2.empty()) return true;
		set<AASequence> best1, best2;
		for (vector<PeptideIdentification>::iterator pep_it = pep1.begin(); pep_it != pep1.end(); ++pep_it)
		{
			if (pep_it->getHits().empty()) continue; // shouldn't be the case
			pep_it->sort();
			best1.insert(pep_it->getHits()[0].getSequence());
		}
		for (vector<PeptideIdentification>::iterator pep_it = pep2.begin(); pep_it != pep2.end(); ++pep_it)
		{
			if (pep_it->getHits().empty()) continue; // shouldn't be the case
			pep_it->sort();
			best2.insert(pep_it->getHits()[0].getSequence());
		}
		return (best1 == best2);
	}

}
