// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
#define OPENMS_ANALYSIS_MAPMATCHING_SIMPLEPAIRPFINDER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>

#if defined OPENMS_DEBUG && ! defined V_SimplePairFinder
#define V_SimplePairFinder(bla) //  std::cout << bla << std::endl;
#else
#define V_SimplePairFinder(bla)
#endif


namespace OpenMS
{

  /**
     @brief This class implements a simple point pair finding algorithm.

     This class implements a feature pair finding algorithm.
     It works on two feature maps, a vector of feature pairs, 
     and a transformation defined for the second feature map (if no
     transformation is given, the pairs are found in the two original maps).

     Policy for copy constructor and assignment: grid_, feature_map_, and
     feature_pairs_ are maintained as pointers and taken shallow copies.  But
     param_ is deep.

  **/
  template < typename MapT = DFeatureMap< 2, DFeature< 2, KernelTraits > > >
  class SimplePairFinder : public BasePairFinder< MapT >
  {
    public:
      typedef DimensionDescription<DimensionDescriptionTagLCMS> DimensionDescriptionType;
      enum DimensionId
      {
        RT = DimensionDescriptionType::RT,
        MZ = DimensionDescriptionType::MZ
    };

      /** @name Type definitions
       */
      //@{
      typedef BasePairFinder< MapT > Base;
      typedef typename Base::TraitsType             TraitsType;

      typedef typename Base::QualityType QualityType;
      typedef typename Base::PositionType           PositionType;
      typedef typename Base::IntensityType          IntensityType;
      typedef typename Base::PointType              PointType;
      typedef typename Base::PointMapType           PointMapType;
      typedef typename Base::FeaturePairType        FeaturePairType;
      typedef typename Base::FeaturePairVectorType  FeaturePairVectorType;
      typedef typename Base::TransformationType     TransformationType;

      using Base::param_;
      using Base::feature_map_;
      using Base::feature_pairs_;
      using Base::transformation_;
      //@}


      ///@name Constructors, destructor and assignment
      //@{
      /// Constructor
      SimplePairFinder()
          : Base(),
          pair_min_quality_(0)
      {}

      /// Copy constructor
      SimplePairFinder(const SimplePairFinder& source)
          : Base(source),
          diff_exponent_(source.diff_exponent_),
          diff_intercept_(source.diff_intercept_),
          pair_min_quality_(source.pair_min_quality_),
          transformed_positions_second_map_(source.transformed_positions_second_map_)
      {}

      ///  Assignment operator
      virtual SimplePairFinder& operator = (SimplePairFinder source)
      {
        param_ = source.param_;
        feature_map_[0] = source.feature_map_[0];
        feature_map_[1] = source.feature_map_[1];
        feature_pairs_ = source.feature_pairs_;
        diff_intercept_[0] = source.diff_intercept_[0];
        diff_intercept_[1] = source.diff_intercept_[1];
        diff_exponent_[0] = source.diff_exponent_[0];
        diff_exponent_[1] = source.diff_exponent_[1];
        pair_min_quality_ = source.pair_min_quality_;
        transformed_positions_second_map_ = source.transformed_positions_second_map_;
        return *this;
      }

      /// Destructor
      virtual ~SimplePairFinder()
      {}
      //@}

      /// returns an instance of this class
      static BasePairFinder<PointMapType>* create()
      {
        return new SimplePairFinder();
      }

      /// returns the name of this module
      static const String getName()
      {
        return "simple";
      }

      template < typename ResultMapType >
      void computeConsensusMap(const PointMapType& first_map, ResultMapType& second_map){}

      /// Estimates the transformation for each grid cell
      virtual void run()
      {
#define V_run(bla) V_SimplePairFinder(bla)
        V_run("@@@ run()");

        parseParam_();

        Size n = feature_map_[1]->size();

        transformed_positions_second_map_.clear();
        transformed_positions_second_map_.resize(n);

        for (Size i = 0; i < n; ++i)
        {
          transformed_positions_second_map_[i] = (*feature_map_[1])[i].getPosition();
        }

        V_run("SimplePairFinder::run(): apply transformation");

        for ( Size dim = 0; dim < 2; ++dim )
        {
          if ( transformation_[dim] )
          {
            for (Size i = 0; i < n; ++i)
            {
              transformation_[dim]->apply( transformed_positions_second_map_[i][dim] );
            }
          }
        }

        V_run("SimplePairFinder::run(): find feature pairs");
#undef V_run

        findFeaturePairs_();
      };

    protected:

      /** @name Data members
       */
      //@{
      /// A parameter for similarity_().
      QualityType diff_exponent_[2];

      /// A parameter for similarity_().
      QualityType diff_intercept_[2];

      /// A parameter for findFeaturePairs_().
      QualityType pair_min_quality_;

      /// The vector of transformed feature positions of the second map
      std::vector<PositionType> transformed_positions_second_map_;

      //@}

      /// Parses the parameters, assigns their values to instance members.
      void parseParam_()
      {
#define V_parseParam_(bla) V_SimplePairFinder(bla)
        V_parseParam_("@@@ parseParam_()");

        {
          std::string param_name_prefix = "similarity:diff_exponent:";
          for ( Size dimension = 0; dimension < 2; ++dimension)
          {
            std::string param_name =
              param_name_prefix + DimensionDescriptionType::dimension_name_short[dimension];
            DataValue data_value = param_.getValue(param_name);
            if ( data_value == DataValue::EMPTY )
            {
              throw Exception::ElementNotFound<std::string>
              (__FILE__,__LINE__,__PRETTY_FUNCTION__,param_name);
            }
            else
            {
              diff_exponent_[dimension] = data_value;
              V_parseParam_(param_name<< ": "<<diff_exponent_[dimension]);
            }
          }
        }

        {
          std::string param_name_prefix = "similarity:diff_intercept:";
          for ( Size dimension = 0; dimension < 2; ++dimension)
          {
            std::string param_name =
              param_name_prefix + DimensionDescriptionType::dimension_name_short[dimension];
            DataValue data_value = param_.getValue(param_name);
            if ( data_value == DataValue::EMPTY )
            {
              throw Exception::ElementNotFound<std::string>
              (__FILE__,__LINE__,__PRETTY_FUNCTION__,param_name);
            }
            else
            {
              diff_intercept_[dimension] = data_value;
              V_parseParam_(param_name<< ": "<<diff_intercept_[dimension]);
            }
          }
        }

        {
          std::string param_name = "similarity:pair_min_quality";
          DataValue data_value = param_.getValue(param_name);
          if ( data_value == DataValue::EMPTY )
          {
            throw Exception::ElementNotFound<std::string>
            (__FILE__,__LINE__,__PRETTY_FUNCTION__,param_name);
          }
          else
          {
            pair_min_quality_ = data_value;
            V_parseParam_(param_name<< ": "<<pair_min_quality_);
          }
        }

#undef V_parseParam_

      } // parseParam_

      /// The actual algorithm for finding feature pairs.
      void findFeaturePairs_()
      {
#define V_findFeaturePairs_(bla) V_SimplePairFinder(bla)
        V_findFeaturePairs_("@@@ findFeaturePairs_()");

        // progress dots
        DataValue const & param_progress_dots = this->param_.getValue("debug:progress_dots");
        int progress_dots
        = param_progress_dots.isEmpty() ? 0 : int(param_progress_dots);
        int number_of_considered_feature_pairs = 0;

        // For each feature in map 0, find his/her best friend in map 1
        std::vector<Size>        best_companion_index_0(feature_map_[0]->size(),Size(-1));
        std::vector<QualityType> best_companion_quality_0(feature_map_[0]->size(),0);
        for ( Size fi0 = 0; fi0 < feature_map_[0]->size(); ++fi0 )
        {
          QualityType best_quality = -std::numeric_limits<QualityType>::max();
          for ( Size fi1 = 0; fi1 < feature_map_[1]->size(); ++ fi1 )
          {
            QualityType quality = similarity_( (*feature_map_[0])[fi0], (*feature_map_[1])[fi1], transformed_positions_second_map_[fi1]);
            if ( quality > best_quality )
            {
              best_quality = quality;
              best_companion_index_0[fi0] = fi1;
            }

            ++number_of_considered_feature_pairs;
            if ( progress_dots &&
                 ! (number_of_considered_feature_pairs % progress_dots)
               )
            {
              std::cout << '-' << std::flush;
            }

          }
          best_companion_quality_0[fi0] = best_quality;
        }

        // For each feature in map 1, find his/her best friend in map 0
        std::vector<Size>        best_companion_index_1(feature_map_[1]->size(),Size(-1));
        std::vector<QualityType> best_companion_quality_1(feature_map_[1]->size(),0);
        for ( Size fi1 = 0; fi1 < feature_map_[1]->size(); ++fi1 )
        {
          QualityType best_quality = -std::numeric_limits<QualityType>::max();
          for ( Size fi0 = 0; fi0 < feature_map_[0]->size(); ++ fi0 )
          {
            QualityType quality = similarity_ ( (*feature_map_[0])[fi0], (*feature_map_[1])[fi1], transformed_positions_second_map_[fi1]);
            if ( quality > best_quality )
            {
              best_quality = quality;
              best_companion_index_1[fi1] = fi0;
            }

            ++number_of_considered_feature_pairs;
            if ( progress_dots &&
                 ! (number_of_considered_feature_pairs % progress_dots)
               )
            {
              std::cout << '+' << std::flush;
            }

          }
          best_companion_quality_1[fi1] = best_quality;
        }

        // And if both like each other, they become a pair.
        // feature_pairs_->clear();
        for ( Size fi0 = 0; fi0 < feature_map_[0]->size(); ++fi0 )
        {
          // fi0 likes someone ...
          if ( best_companion_quality_0[fi0] > pair_min_quality_ )
          {
            // ... who likes him too ...
            Size best_companion_of_fi0 = best_companion_index_0[fi0];
            if ( best_companion_index_1[best_companion_of_fi0] == fi0 &&
                 best_companion_quality_1[best_companion_of_fi0] > pair_min_quality_
               )
            {
              feature_pairs_->push_back
              ( FeaturePairType ( (*feature_map_[0])[fi0],
                                  (*feature_map_[1])[best_companion_of_fi0],
                                  best_companion_quality_0[fi0] + best_companion_quality_1[best_companion_of_fi0]
                                )
              );
            }
          }
        }

#undef V_findFeaturePairs_

      } // findFeaturePairs_

      /**@brief Compute the similarity for a pair of features; larger quality
        values are better.

        The returned value should express our confidence that one feature might
        possibly be matched to the other.

        The details here are kind of alchemy ...
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

        for ( Size dimension = 0; dimension < 2; ++dimension )
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

        return intensity_ratio / position_difference[RT] / position_difference[MZ];
      }

  }
  ; // SimplePairFinder

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_SimplePairFinder_H
