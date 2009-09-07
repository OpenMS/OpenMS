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

namespace OpenMS
{

  StablePairFinder::StablePairFinder() :
    Base()
  {
    //set the name for DefaultParamHandler error messages
    Base::setName(getProductName());

    defaults_.setValue("diff_exponent:RT", 1.0,
      "RT differences are raised to this power", StringList::create("advanced"));
    defaults_.setMinFloat("diff_exponent:RT",0.);

    defaults_.setValue("diff_exponent:MZ", 2.0,
      "MZ differences are raised to this power", StringList::create("advanced"));
    defaults_.setMinFloat("diff_exponent:MZ",0.);

    defaults_.setSectionDescription("diff_exponent",
      "Absolute position differences are raised to this power. "
        "E.g. 1 for 'linear' distance, 2 for 'quadratic' distance");

    defaults_.setValue("intensity_exponent", 0.5,
      "Intensity ratios are raised to this power.  "
        "If set to 0, intensities are not considered.", StringList::create(
        "advanced"));
    defaults_.setMinFloat("intensity_exponent",0.);

    defaults_.setValue("max_pair_distance:RT", 100.0,
      "Maximal allowed distance in RT for a pair, when MZ is equal");
    defaults_.setMinFloat("max_pair_distance:RT",0.);

    defaults_.setValue("max_pair_distance:MZ", 0.3,
      "Maximal allowed distance in MZ for a pair, when RT is equal");
    defaults_.setMinFloat("max_pair_distance:MZ",0.);

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

    Base::defaultsToParam_();
  }

  void
  StablePairFinder::updateMembers_()
  {
    V_("@@@ StablePairFinder::updateMembers_()");

    diff_exponent_[RT] = (DoubleReal) param_.getValue("diff_exponent:RT");
    if ( diff_exponent_[RT] <= 0 )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__,
        "exponent for RT distance must be > 0");
    }
    VV_(diff_exponent_[RT]);

    diff_exponent_[MZ] = (DoubleReal) param_.getValue("diff_exponent:MZ");
    if ( diff_exponent_[MZ] <= 0 )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__,
        "exponent for MZ distance must be > 0");
    }
    VV_(diff_exponent_[MZ]);

    intensity_exponent_ = (DoubleReal) param_.getValue("intensity_exponent");
    if ( intensity_exponent_ < 0 )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__,
        "exponent for intensity ratio must be >= 0");
    }
    VV_(intensity_exponent_);

    max_pair_distance_[MZ] = (DoubleReal) param_.getValue(
      "max_pair_distance:MZ");
    if ( max_pair_distance_[MZ] < 0 )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__,
        "max pair distance for MZ must be >= 0");
    }
    VV_(max_pair_distance_[MZ]);
    max_pair_distance_reciprocal_[MZ] = 1. / max_pair_distance_[MZ];
    VV_(max_pair_distance_reciprocal_[MZ]);

    max_pair_distance_[RT] = (DoubleReal) param_.getValue(
      "max_pair_distance:RT");
    if ( max_pair_distance_[RT] < 0 )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__,
        "max pair distance for RT must be >= 0");
    }
    VV_(max_pair_distance_[RT]);
    max_pair_distance_reciprocal_[RT] = 1. / max_pair_distance_[RT];
    VV_(max_pair_distance_reciprocal_[RT]);

    different_charge_penalty_ = (DoubleReal) param_.getValue("different_charge_penalty");
    if ( different_charge_penalty_ < 1. )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__,
        "different_charge_penalty must be >= 1");
    }
    VV_(different_charge_penalty_);

    second_nearest_gap_ = (DoubleReal) param_.getValue("second_nearest_gap");
    if ( second_nearest_gap_ <= 1 )
    {
      throw Exception::InvalidParameter(__FILE__,__LINE__, __PRETTY_FUNCTION__,
        "second_nearest_gap must be >= 1");
    }
    VV_(second_nearest_gap_);

    return;
  }

  DoubleReal
  StablePairFinder::distance_( ConsensusFeature const & left,
                               ConsensusFeature const & right ) const
  {
    DPosition<2> position_difference = left.getPosition() - right.getPosition();
    for ( UInt dimension = 0; dimension < 2; ++dimension )
    {
      if ( position_difference[dimension] < 0 )
      {
        position_difference[dimension] = -position_difference[dimension];
      }
      position_difference[dimension]
          *= max_pair_distance_reciprocal_[dimension];
      position_difference[dimension] = pow(position_difference[dimension],
        diff_exponent_[dimension]);
    }
    DoubleReal result = position_difference[RT] + position_difference[MZ];
    if ( intensity_exponent_ != 0 )
    {
      DoubleReal right_intensity(right.getIntensity());
      DoubleReal left_intensity(left.getIntensity());
      if ( right_intensity == 0 || left_intensity == 0 )
      {
        return std::numeric_limits<DoubleReal>::max();
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

    // For each element in map 0, find its best friend in map 1
    std::vector<UInt> best_companion_index_0(input_maps[0].size(), UInt(-1));
    std::vector<DoubleReal> best_companion_distance_0(input_maps[0].size(), 0);
    std::vector<DoubleReal> second_best_companion_distance_0(
      input_maps[0].size(), 0);
    for ( UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0 )
    {
      DoubleReal best_distance = std::numeric_limits<DoubleReal>::max();
      DoubleReal second_best_distance = std::numeric_limits<DoubleReal>::max();
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
    std::vector<UInt> best_companion_index_1(input_maps[1].size(), UInt(-1));
    std::vector<DoubleReal> best_companion_distance_1(input_maps[1].size(), 0);
    std::vector<DoubleReal> second_best_companion_distance_1(
      input_maps[1].size(), 0);
    for ( UInt fi1 = 0; fi1 < input_maps[1].size(); ++fi1 )
    {
      DoubleReal best_distance = std::numeric_limits<DoubleReal>::max();
      DoubleReal second_best_distance = std::numeric_limits<DoubleReal>::max();
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
      // fi0 likes someone ...
      if ( best_companion_distance_0[fi0] < 1. && best_companion_distance_0[fi0]
          * second_nearest_gap_ <= second_best_companion_distance_0[fi0] )
      {
        // ... who likes him too ...
        UInt best_companion_of_fi0 = best_companion_index_0[fi0];
        if ( best_companion_index_1[best_companion_of_fi0] == fi0
            && best_companion_distance_1[best_companion_of_fi0]
                * second_nearest_gap_
                <= second_best_companion_distance_1[best_companion_of_fi0] )
        {
          result_map.push_back(ConsensusFeature());
          ConsensusFeature & f = result_map.back();

          f.insert(input_maps[0][fi0]);
          f.getPeptideIdentifications().
            insert( f.getPeptideIdentifications().end(),
              input_maps[0][fi0].getPeptideIdentifications().begin(),
              input_maps[0][fi0].getPeptideIdentifications().end()
            );

          f.insert(input_maps[1][best_companion_of_fi0]);
          f.getPeptideIdentifications().
            insert( f.getPeptideIdentifications().end(),
              input_maps[1][best_companion_of_fi0].getPeptideIdentifications().begin(),
              input_maps[1][best_companion_of_fi0].getPeptideIdentifications().end()
            );

          f.computeConsensus();
          DoubleReal quality = 1. - best_companion_distance_0[fi0];
          DoubleReal quality0 = 1. - best_companion_distance_0[fi0]
              * second_nearest_gap_ / second_best_companion_distance_0[fi0];
          DoubleReal quality1 = 1.
              - best_companion_distance_1[best_companion_of_fi0]
                  * second_nearest_gap_
                  / second_best_companion_distance_1[best_companion_of_fi0];
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

}
