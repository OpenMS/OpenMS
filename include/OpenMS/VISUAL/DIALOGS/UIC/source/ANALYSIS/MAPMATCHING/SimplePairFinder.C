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

#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>

namespace OpenMS
{
	SimplePairFinder::SimplePairFinder()
		: Base()
  {
  	//set the name for DefaultParamHandler error messages
		setName(getProductName());

		defaults_.setValue("similarity:diff_intercept:RT",1.0,"This parameter controls the asymptotic decay rate for large differences (for more details see the similarity measurement).",StringList::create("advanced"));
    defaults_.setValue("similarity:diff_intercept:MZ",0.1,"This parameter controls the asymptotic decay rate for large differences (for more details see the similarity measurement).",StringList::create("advanced"));
    defaults_.setValue("similarity:diff_exponent:RT",2.0,"This parameter is important for small differences (for more details see the similarity measurement).",StringList::create("advanced"));
    defaults_.setValue("similarity:diff_exponent:MZ",1.0,"This parameter is important for small differences (for more details see the similarity measurement).",StringList::create("advanced"));
    defaults_.setValue("similarity:pair_min_quality",0.01,"Minimum required pair quality.",StringList::create("advanced"));

    Base::defaultsToParam_();
  }

	void SimplePairFinder::run(const std::vector<ConsensusMap>& input_maps, ConsensusMap& result_map)
  {
  	if (input_maps.size()!=2) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"exactly two input maps required");
		checkIds_(input_maps);

    // progress dots
    Int progress_dots = 0;
		if (this->param_.exists("debug::progress_dots"))
		{
    	progress_dots = (Int)this->param_.getValue("debug:progress_dots");
	 	}
		Int number_of_considered_element_pairs = 0;

    // For each element in map 0, find its best friend in map 1
    std::vector<UInt>        best_companion_index_0(input_maps[0].size(),UInt(-1));
    std::vector<DoubleReal>  best_companion_quality_0(input_maps[0].size(),0);
    for ( UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0 )
		{
			DoubleReal best_quality = -std::numeric_limits<DoubleReal>::max();
			for ( UInt fi1 = 0; fi1 < input_maps[1].size(); ++ fi1 )
			{
				DoubleReal quality = similarity_( input_maps[0][fi0], input_maps[1][fi1]);
				if ( quality > best_quality )
				{
					best_quality = quality;
					best_companion_index_0[fi0] = fi1;
				}

				++number_of_considered_element_pairs;
				if ( progress_dots && ! (number_of_considered_element_pairs % progress_dots) )
				{
					std::cout << '-' << std::flush;
				}

			}
			best_companion_quality_0[fi0] = best_quality;
    }

		// For each element in map 1, find its best friend in map 0
		std::vector<UInt>       best_companion_index_1(input_maps[1].size(),UInt(-1));
    std::vector<DoubleReal> best_companion_quality_1(input_maps[1].size(),0);
    for ( UInt fi1 = 0; fi1 < input_maps[1].size(); ++fi1 )
		{
			DoubleReal best_quality = -std::numeric_limits<DoubleReal>::max();
			for ( UInt fi0 = 0; fi0 < input_maps[0].size(); ++ fi0 )
			{
				DoubleReal quality = similarity_( input_maps[0][fi0], input_maps[1][fi1]);
				if ( quality > best_quality )
				{
					best_quality = quality;
					best_companion_index_1[fi1] = fi0;
				}

				++number_of_considered_element_pairs;
				if ( progress_dots &&
						 ! (number_of_considered_element_pairs % progress_dots)
					 )
				{
					std::cout << '+' << std::flush;
				}

			}
			best_companion_quality_1[fi1] = best_quality;
    }

    // And if both like each other, they become a pair.
    // element_pairs_->clear();
		for ( UInt fi0 = 0; fi0 < input_maps[0].size(); ++fi0 )
    {
			// fi0 likes someone ...
			if ( best_companion_quality_0[fi0] > pair_min_quality_ )
			{
				// ... who likes him too ...
				UInt best_companion_of_fi0 = best_companion_index_0[fi0];
				if ( best_companion_index_1[best_companion_of_fi0] == fi0 &&
						 best_companion_quality_1[best_companion_of_fi0] > pair_min_quality_
					 )
				{
					ConsensusFeature f;
					f.insert( input_maps[0][fi0] );
					f.insert( input_maps[1][best_companion_of_fi0] );
					f.computeConsensus();
					f.setQuality(best_companion_quality_0[fi0] + best_companion_quality_1[best_companion_of_fi0]);
					result_map.push_back(f);
				}
			}
    }
		return;
  }

  void SimplePairFinder::updateMembers_()
  {
    diff_intercept_[Peak2D::RT] = (DoubleReal)param_.getValue("similarity:diff_intercept:RT");
		if ( diff_intercept_[Peak2D::RT] <= 0 )
		{
			throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"intercept for RT must be > 0");
		}

    diff_intercept_[Peak2D::MZ] = (DoubleReal)param_.getValue("similarity:diff_intercept:MZ");
		if ( diff_intercept_[Peak2D::MZ] <= 0 )
		{
			throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__,"intercept for MZ must be > 0");
		}

    diff_exponent_[Peak2D::RT] = (DoubleReal)param_.getValue("similarity:diff_exponent:RT");
    diff_exponent_[Peak2D::MZ] = (DoubleReal)param_.getValue("similarity:diff_exponent:MZ");
    pair_min_quality_ = (DoubleReal)param_.getValue("similarity:pair_min_quality");
  }

  DoubleReal SimplePairFinder::similarity_ ( ConsensusFeature const & left, ConsensusFeature const & right) const
  {
    DoubleReal right_intensity(right.getIntensity());
    if ( right_intensity == 0 ) return 0;
    DoubleReal intensity_ratio = left.getIntensity() / right_intensity;
    if ( intensity_ratio > 1. ) intensity_ratio = 1. / intensity_ratio;

    // if the right map is the transformed map, take the transformed right position
    DPosition<2> position_difference = left.getPosition() - right.getPosition();

    for ( UInt dimension = 0; dimension < 2; ++dimension )
    {
			// the formula is explained in class doc
			if ( position_difference[dimension] < 0 )
			{
				position_difference[dimension] = -position_difference[dimension];
			}
			position_difference[dimension] *= diff_intercept_[dimension];
			position_difference[dimension] += 1.0;
			position_difference[dimension] = pow(position_difference[dimension],diff_exponent_[dimension]);
    }

    return intensity_ratio / position_difference[Peak2D::RT] / position_difference[Peak2D::MZ];
  }

}
