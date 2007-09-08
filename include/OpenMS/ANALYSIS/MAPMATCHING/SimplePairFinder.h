// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

     @note This pair finder does not offer a method to compute consensus elements given
     two element maps!

		 @ref SimplePairFinder_Parameters are explained on a separate page.
  */
  template < typename MapT = FeatureMap< > >
  class SimplePairFinder : public BasePairFinder< MapT >
  {
	 public:
    /** Symbolic names for indices of element maps etc.
				This should make things more understandable and maintainable.
		*/
    enum Maps
			{
				MODEL = 0,
				SCENE = 1
			};

    typedef BasePairFinder< MapT > Base;
    //typedef typename Base::TraitsType             TraitsType;

    typedef typename Base::QualityType QualityType;
    typedef typename Base::PositionType           PositionType;
    typedef typename Base::IntensityType          IntensityType;
    typedef typename Base::PointType              PointType;
    typedef typename Base::PointMapType           PointMapType;
    typedef typename Base::ElementPairType        ElementPairType;
    typedef typename Base::ElementPairVectorType  ElementPairVectorType;
    typedef typename Base::TransformationType     TransformationType;

    using Base::param_;
    using Base::defaults_;
    using Base::setName;
    using Base::element_map_;
    using Base::element_pairs_;
    using Base::transformation_;

    /// Constructor
    SimplePairFinder()
			: Base()
    {
    	//set the name for DefaultParamHandler error messages
			setName(getProductName());
			
			defaults_.setValue("similarity:diff_intercept:RT",1,"This parameter controls the asymptotic decay rate for large differences (for more details see the similarity measurement).");
      defaults_.setValue("similarity:diff_intercept:MZ",0.1,"This parameter controls the asymptotic decay rate for large differences (for more details see the similarity measurement).");
      defaults_.setValue("similarity:diff_exponent:RT",2,"This parameter is important for small differences (for more details see the similarity measurement).");
      defaults_.setValue("similarity:diff_exponent:MZ",1,"This parameter is important for small differences (for more details see the similarity measurement).");
      defaults_.setValue("similarity:pair_min_quality",0.01,"Minimum required pair quality.");

      Base::defaultsToParam_();
    }

    /// Destructor
    virtual ~SimplePairFinder()
		{
		}

    /// returns an instance of this class
    static BasePairFinder<PointMapType>* create()
    {
      return new SimplePairFinder();
    }

    /// returns the name of this module
    static const String getProductName()
    {
      return "SimplePairFinder";
    }

    /// Get diff exponent. See @sa similarity_().
    double getDiffExponent(UInt dim)
    {
      return diff_exponent_[dim];
    }

    /// Set diff exponent. See @sa similarity_().
    void setDiffExponent(UInt dim, DoubleReal exponent)
    {
      diff_exponent_[dim] = exponent;
      String param_name_prefix = "similarity:diff_exponent:";
      String param_name = param_name_prefix + RawDataPoint2D::shortDimensionName(dim);
      param_.setValue(param_name, exponent);
    }

    /// Get diff intercept. See @sa similarity_().
    double getDiffIntercept(UInt dim)
    {
      return diff_intercept_[dim];
    }

    /// Set diff intercept. See @sa similarity_().
    void setDiffIntercept(UInt dim, DoubleReal intercept)
    {
      diff_intercept_[dim] = intercept;
      param_.setValue(String("similarity:diff_intercept:") + RawDataPoint2D::shortDimensionName(dim), intercept);
    }

    /// Get pair min quality
    double getPairMinQuality()
    {
      return pair_min_quality_;
    }

    /// Set pair min quality
    void setPairMinQuality(DoubleReal quality)
    {
      pair_min_quality_ = quality;
      param_.setValue("similarity:pair_min_quality", quality);
    }

    /** @brief Find pairs of elements in both maps.

    For each feature, we find the nearest neighbor in the other map according to @sa similarity_().
    If two features point at each other, they become a pair.
    */
    virtual void findElementPairs()
    {
#define V_findElementPairs(bla) V_SimplePairFinder(bla)
      UInt n = element_map_[SCENE]->size();

      transformed_positions_second_map_.clear();
      transformed_positions_second_map_.resize(n);

      for (UInt i = 0; i < n; ++i)
      {
				transformed_positions_second_map_[i] = (*element_map_[SCENE])[i].getPosition();
      }

      V_findElementPairs("SimplePairFinder::run(): apply transformation");

      for ( UInt dim = 0; dim < 2; ++dim )
      {
				for (UInt i = 0; i < n; ++i)
				{
					transformation_[dim].apply( transformed_positions_second_map_[i][dim] );
				}
      }

      V_findElementPairs("SimplePairFinder::run(): find element pairs" << pair_min_quality_);

      // progress dots
      Int progress_dots = 0;
			if (this->param_.exists("debug::progress_dots"))
			{
      	progress_dots = (Int)this->param_.getValue("debug:progress_dots");
		 	}
			Int number_of_considered_element_pairs = 0;

      // For each element in map 0, find his/her best friend in map 1
      std::vector<UInt>        best_companion_index_0(element_map_[MODEL]->size(),UInt(-1));
      std::vector<QualityType> best_companion_quality_0(element_map_[MODEL]->size(),0);
      for ( UInt fi0 = 0; fi0 < element_map_[MODEL]->size(); ++fi0 )
      {
				QualityType best_quality = -std::numeric_limits<QualityType>::max();
				for ( UInt fi1 = 0; fi1 < element_map_[SCENE]->size(); ++ fi1 )
				{
					QualityType quality = similarity_( (*element_map_[MODEL])[fi0], (*element_map_[SCENE])[fi1], transformed_positions_second_map_[fi1]);
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
      std::vector<UInt>        best_companion_index_1(element_map_[SCENE]->size(),UInt(-1));
      std::vector<QualityType> best_companion_quality_1(element_map_[SCENE]->size(),0);
      for ( UInt fi1 = 0; fi1 < element_map_[SCENE]->size(); ++fi1 )
      {
				QualityType best_quality = -std::numeric_limits<QualityType>::max();
				for ( UInt fi0 = 0; fi0 < element_map_[MODEL]->size(); ++ fi0 )
				{
					QualityType quality = similarity_( (*element_map_[MODEL])[fi0], (*element_map_[SCENE])[fi1], transformed_positions_second_map_[fi1]);
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
      for ( UInt fi0 = 0; fi0 < element_map_[MODEL]->size(); ++fi0 )
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
						element_pairs_->push_back( ElementPairType ( (*element_map_[SCENE])[best_companion_of_fi0],
																												 (*element_map_[MODEL])[fi0],
																												 best_companion_quality_0[fi0] + best_companion_quality_1[best_companion_of_fi0]
																											 ));
					}
				}
      }

#undef V_findElementPairs

    } // findElementPairs


	 protected:
    virtual void updateMembers_()
    {
      diff_intercept_[RawDataPoint2D::RT] = (QualityType)param_.getValue("similarity:diff_intercept:RT");
      diff_intercept_[RawDataPoint2D::MZ] = (QualityType)param_.getValue("similarity:diff_intercept:MZ");
      diff_exponent_[RawDataPoint2D::RT] = (QualityType)param_.getValue("similarity:diff_exponent:RT");
      diff_exponent_[RawDataPoint2D::MZ] = (QualityType)param_.getValue("similarity:diff_exponent:MZ");
      pair_min_quality_ = (QualityType)param_.getValue("similarity:pair_min_quality");
    }

    /// A parameter for #similarity_().
    QualityType diff_exponent_[2];

    /// A parameter for #similarity_().
    QualityType diff_intercept_[2];

    /// A parameter for findElementPairs_().
    QualityType pair_min_quality_;

    /// The vector of transformed element positions of the second map
    std::vector<PositionType> transformed_positions_second_map_;

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
    QualityType similarity_ ( PointType const & left, PointType const & right, const PositionType& new_position) const
    {
      QualityType right_intensity(right.getIntensity());
      if ( right_intensity == 0 )
				return 0;
      QualityType intensity_ratio = left.getIntensity() / right_intensity;
      if ( intensity_ratio > 1. )
				intensity_ratio = 1. / intensity_ratio;

      // if the right map is the transformed map, take the transformed right position
      PositionType position_difference = left.getPosition() - new_position;

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
