// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_SIMPLEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_SIMPLEPAIRFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>

#define V_SimplePairFinder(bla) // std::cout << bla << std::endl;

namespace OpenMS
{

  /**
		@brief This class implements a simple point pair finding algorithm.
		
		It offers a method to determine element pairs in two element maps,
		given two point maps and a transformation defined for the second element map (if no
		transformation is given, the pairs are found in the two original maps).
		
		@ref SimplePairFinder_Parameters are explained on a separate page.

		@ingroup FeatureGrouping
  */
  class SimplePairFinder
  	: public BasePairFinder
  {
	 public:
	 	///Base class
    typedef BasePairFinder Base;

    /// Constructor
    SimplePairFinder()
			: Base()
    {
    	//set the name for DefaultParamHandler error messages
			setName(getProductName());
			
			defaults_.setValue("similarity:diff_intercept:RT",1.0,"This parameter controls the asymptotic decay rate for large differences (for more details see the similarity measurement).",true);
      defaults_.setValue("similarity:diff_intercept:MZ",0.1,"This parameter controls the asymptotic decay rate for large differences (for more details see the similarity measurement).",true);
      defaults_.setValue("similarity:diff_exponent:RT",2.0,"This parameter is important for small differences (for more details see the similarity measurement).",true);
      defaults_.setValue("similarity:diff_exponent:MZ",1.0,"This parameter is important for small differences (for more details see the similarity measurement).",true);
      defaults_.setValue("similarity:pair_min_quality",0.01,"Minimum required pair quality.",true);

      Base::defaultsToParam_();
    }

    /// Destructor
    virtual ~SimplePairFinder()
		{
		}

    /// returns an instance of this class
    static BasePairFinder* create()
    {
      return new SimplePairFinder();
    }

    /// returns the name of this module
    static const String getProductName()
    {
      return "simple";
    }

    /** 
    	@brief Find pairs of elements in both maps.

	    For each feature, we find the nearest neighbor in the other map according to @sa similarity_().
	    If two features point at each other, they become a pair.
    */
		virtual void run(ConsensusMap& result_map)
    {
      UInt n = scene_map_->size();

      transformed_positions_second_map_.clear();
      transformed_positions_second_map_.resize(n);

      for (UInt i = 0; i < n; ++i)
      {
				transformed_positions_second_map_[i] = (*scene_map_)[i].getPosition();
      }
      
      // progress dots
      Int progress_dots = 0;
			if (this->param_.exists("debug::progress_dots"))
			{
      	progress_dots = (Int)this->param_.getValue("debug:progress_dots");
		 	}
			Int number_of_considered_element_pairs = 0;

      // For each element in map 0, find his/her best friend in map 1
      std::vector<UInt>        best_companion_index_0(getModelMap().size(),UInt(-1));
      std::vector<DoubleReal> best_companion_quality_0(getModelMap().size(),0);
      for ( UInt fi0 = 0; fi0 < getModelMap().size(); ++fi0 )
			{
				DoubleReal best_quality = -std::numeric_limits<DoubleReal>::max();
				for ( UInt fi1 = 0; fi1 < getSceneMap().size(); ++ fi1 )
				{
					DoubleReal quality = similarity_( getModelMap()[fi0], getSceneMap()[fi1], transformed_positions_second_map_[fi1]);
					if ( quality > best_quality )
					{
						best_quality = quality;
						best_companion_index_0[fi0] = fi1;
					}

					++number_of_considered_element_pairs;
					if ( progress_dots &&
							 ! (number_of_considered_element_pairs % progress_dots)
						 )
					{
						std::cout << '-' << std::flush;
					}

				}
				best_companion_quality_0[fi0] = best_quality;
      }

			// For each element in map 1, find his/her best friend in map 0
			std::vector<UInt>        best_companion_index_1(getSceneMap().size(),UInt(-1));
      std::vector<DoubleReal> best_companion_quality_1(getSceneMap().size(),0);
      for ( UInt fi1 = 0; fi1 < getSceneMap().size(); ++fi1 )
			{
				DoubleReal best_quality = -std::numeric_limits<DoubleReal>::max();
				for ( UInt fi0 = 0; fi0 < getModelMap().size(); ++ fi0 )
				{
					DoubleReal quality = similarity_( getModelMap()[fi0], getSceneMap()[fi1], transformed_positions_second_map_[fi1]);
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
			for ( UInt fi0 = 0; fi0 < getModelMap().size(); ++fi0 )
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
						if ( model_index_ == -1)
						{
							f.insert( getModelMap()[fi0] );
						}
						else
						{
							f.insert( model_index_, fi0, getModelMap()[fi0] );
						}
						if ( scene_index_ == -1)
						{
							f.insert( getSceneMap()[best_companion_of_fi0] );
						}
						else
						{
							f.insert( scene_index_, best_companion_of_fi0, getSceneMap()[best_companion_of_fi0] );
						}
						f.computeConsensus();
						f.setQuality(best_companion_quality_0[fi0] + best_companion_quality_1[best_companion_of_fi0]);
						result_map.push_back(f);
					}
				}
      }
    }


	 protected:
    virtual void updateMembers_()
    {
      diff_intercept_[RawDataPoint2D::RT] = (DoubleReal)param_.getValue("similarity:diff_intercept:RT"); // TODO: use internal_mz_scaling instead
      diff_intercept_[RawDataPoint2D::MZ] = (DoubleReal)param_.getValue("similarity:diff_intercept:MZ"); // TODO: use internal_mz_scaling instead
      diff_exponent_[RawDataPoint2D::RT] = (DoubleReal)param_.getValue("similarity:diff_exponent:RT");
      diff_exponent_[RawDataPoint2D::MZ] = (DoubleReal)param_.getValue("similarity:diff_exponent:MZ");
      pair_min_quality_ = (DoubleReal)param_.getValue("similarity:pair_min_quality");
    }

    /// A parameter for similarity_().
    DoubleReal diff_exponent_[2];

    /// A parameter for similarity_().
		 // TODO: use internal_mz_scaling instead
    DoubleReal diff_intercept_[2];

    /// Minimal pair quality
    DoubleReal pair_min_quality_;

    /// The vector of transformed element positions of the second map
    std::vector<DPosition<2> > transformed_positions_second_map_;

    /**@brief Compute the similarity for a pair of elements; larger quality
		values are better.

		The returned value should express our confidence that one element might
		possibly be matched to the other.

		The similarity is computed as follows.
		- For each dimension:
		<ul>
		<li> Take the absolute difference of the coordinates.
		<li> Add #diff_intercept_ to it.
		<li> Raise the sum to power of #diff_exponent_.
		</ul>
		- Multiply these numbers for both dimensions.
		- Take the reciprocal value of the result.
		.

		The parameter #diff_exponent_ controls the asymptotic decay rate for large
		differences.  The parameter #diff_intercept_ is important for small
		differences.

    */
		// TODO do this with FeatureHandle instead, dont hand over new_position separately
    DoubleReal similarity_ ( ConsensusFeature const & left, ConsensusFeature const & right, const DPosition<2>& new_position) const
    {
      DoubleReal right_intensity(right.getIntensity());
      if ( right_intensity == 0 )
				return 0;
      DoubleReal intensity_ratio = left.getIntensity() / right_intensity;
      if ( intensity_ratio > 1. )
				intensity_ratio = 1. / intensity_ratio;

      // if the right map is the transformed map, take the transformed right position
      DPosition<2> position_difference = left.getPosition() - new_position;

      for ( UInt dimension = 0; dimension < 2; ++dimension )
      {
				// Take the absolute value
				if ( position_difference[dimension] < 0 )
				{
					position_difference[dimension] = -position_difference[dimension];
				}
				// Raise the difference to a (potentially fractional) power
				position_difference[dimension] =
					pow(position_difference[dimension],diff_exponent_[dimension]);
				// Add an absolute number
				position_difference[dimension] += diff_intercept_[dimension];
      }

      return intensity_ratio / position_difference[RawDataPoint2D::RT] / position_difference[RawDataPoint2D::MZ];
    }

  }
  ; // SimplePairFinder

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_SimplePairFinder_H
