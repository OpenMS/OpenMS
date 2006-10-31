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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_GEOMHASHSHIFTSUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_GEOMHASHSHIFTSUPERIMPOSER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DLinearMapping.h>


#include <fstream>

#if defined OPENMS_DEBUG && ! defined V_GeomHashShiftSuperimposer
#define V_GeomHashShiftSuperimposer(bla) // std::cout << bla << std::endl;
#else
#define V_GeomHashShiftSuperimposer(bla)
#endif

namespace OpenMS
{

  /**
     @brief Superimposer that uses geometric hashing to find a good shift

     While the feature positions can have D dimensions, only the first two are
     used to find a shift.
  **/
  template < typename MapT = DFeatureMap<2> >
  class GeomHashShiftSuperimposer
        : public BaseSuperimposer< MapT >
  {
    public:

      /** @name Type definitions
       */
      //@{
      /// Defines the coordinates of peaks / features.
      typedef DimensionDescription<DimensionDescriptionTagLCMS> DimensionDescriptionType;
      /// Defines the coordinates of peaks / features.
      enum DimensionId
      {
        RT = DimensionDescriptionType::RT,
        MZ = DimensionDescriptionType::MZ
    };

      /** Symbolic names for indices of feature maps etc.
          This should make things more understandable and maintainable.
           */
      enum Maps
      {
        MODEL = 0,
        SCENE = 1
    };

      typedef BaseSuperimposer< MapT > Base;
      //@}

    protected:
      // We need this to make the intensity bounding box use the intensity type
      // instead of the coordinate type.
    struct IntensityBoundingBoxTraits : Base::TraitsType
      {
        typedef typename Base::TraitsType::IntensityType CoordinateType;
      };

    public:

      /** @brief Nested class to represent a shift.

      The shift itself is stored as a DPosition.  Also provided is a
      quality value, with an acompanying comparator.
      */
      class Shift
      {
        public:

          Shift()
              : position_(0),
              quality_(0)
          {}

          Shift(Shift const & source)
              : position_(source.position_),
              quality_(source.quality_)
          {}

          Shift & operator= (Shift const & source)
          {
            position_ = source.position_;
            quality_ = source.quality_;
            return *this;
          }

          ~Shift()
          {}

          typedef typename MapT::TraitsType TraitsType;
          typedef DPosition<2,TraitsType> PositionType;
          typedef typename TraitsType::QualityType QualityType;

          /// Non-mutable access to the data point position (multidimensional)
          const PositionType& getPosition() const
          {
            return position_;
          }

          /// Mutable access to the data point position (multidimensional)
          PositionType& getPosition()
          {
            return position_;
          }

          /// Mutable access to the data point position (multidimensional)
          void setPosition(const PositionType& position)
          {
            position_ = position;
          }

          /// Non-mutable access to the quality
          const QualityType& getQuality() const
          {
            return quality_;
          }

          /// Mutable access to the quality
          QualityType& getQuality()
          {
            return quality_;
          }

          /// Mutable access to the quality
          void setQuality(const QualityType& quality)
          {
            quality_ = quality;
          }

          /*

          not used currently

          /// Compare by getQuality()
          struct QualityLess : std::binary_function < Shift, Shift, bool >
          {
          bool operator () ( Shift const & left, Shift const & right ) const {return ( left.getQuality() < right.getQuality() );}
          bool operator () ( Shift const & left, QualityType const & right ) const {return ( left.getQuality() < right );}
          bool operator () ( QualityType const & left, Shift const & right ) const {return ( left< right.getQuality() );}
          bool operator () ( QualityType const & left, QualityType const & right ) const {return ( left < right );}
          };
          */

        protected:
          PositionType position_;
          QualityType quality_;
      };

      /** @name Type definitions
       */
      //@{
      typedef typename Base::TraitsType TraitsType;
      typedef typename Base::QualityType QualityType;
      typedef typename Base::PositionType PositionType;
      typedef typename Base::IntensityType IntensityType;
      typedef typename Base::PointType PointType;
      typedef typename Base::PointMapType PointMapType;
      typedef DBoundingBox<2,TraitsType>  PositionBoundingBoxType;
      typedef DBoundingBox<1,IntensityBoundingBoxTraits> IntensityBoundingBoxType;
      typedef std::vector <Size> FeatureBucketType;
      typedef Matrix < FeatureBucketType > FeatureBucketMatrixType;
      typedef Shift ShiftType;
      typedef Matrix < typename ShiftType::QualityType > ShiftQualityMatrixType;
      typedef Matrix < ShiftType > ShiftMatrixType;
      typedef DLinearMapping< 1, TraitsType > FinalShiftType;
      //@}
      using Base::setParam;
      using Base::getParam;
      using Base::setFeatureMap;
      using Base::getFeatureMap;
      using Base::final_transformation_;

      //------------------------------------------------------------

      ///@name Constructors, destructor and assignment
      //@{
      /// Constructor
      GeomHashShiftSuperimposer()
          : Base()
      {
        final_transformation_[RT] = new FinalShiftType();
        final_transformation_[MZ] = new FinalShiftType();
      }

      /// Copy constructor
      GeomHashShiftSuperimposer(const GeomHashShiftSuperimposer& source)
          : Base(source)
      {
        final_transformation_[RT] = new FinalShiftType(source.final_transformation_[RT]);
        final_transformation_[MZ] = new FinalShiftType(source.final_transformation_[MZ]);
      }

      ///  Assignment operator
      GeomHashShiftSuperimposer& operator = (const GeomHashShiftSuperimposer& source)
      {
        Base::operator=(source);
        for ( Size dim = 0; dim < 2; ++dim )
        {
          if (final_transformation_[dim] != 0)
          {
            final_transformation_[dim] = new FinalShiftType(source.final_transformation_[dim]);
          }
          else
          {
            *final_transformation_[dim] = *(source.final_transformation_[dim]);
          }
        }
        return *this;
      }

      /// Destructor
      virtual ~GeomHashShiftSuperimposer()
      {
        V_GeomHashShiftSuperimposer("~GeomHashShift");
      }
      //@}

      //------------------------------------------------------------

      /** @name Accesssor methods
       */
      /// Set
      virtual void setShift(Size dim, const FinalShiftType& shift)
      {
        if (final_transformation_[dim] != 0)
        {
          final_transformation_[dim] = new FinalShiftType(shift);
        }
      }

      /// Get shift in dimension dim
      virtual const FinalShiftType& getShift(Size dim) const
      {
        FinalShiftType* f;
        if ((f = dynamic_cast<FinalShiftType*>(final_transformation_[dim])))
        {
          return *f;
        }
        else
        {
          // do what ????
          exit (99);
        }
      }

      //----------------------------------------------------------------------

      /// estimates the transformation for each grid cell
      virtual void run()
      {
        /// clear the member
        feature_bucket_[RT].clear();
        feature_bucket_[MZ].clear();
        shift_bucket_.clear();
        shift_matrix_.clear();

        for ( Size dim = 0; dim < 2; ++dim )
        {
          if (final_transformation_[dim] != 0)
          {
            delete final_transformation_[dim];
            final_transformation_[dim] = 0;
          }
        }

        if ( !this->feature_map_[MODEL]->empty() && !this->feature_map_[SCENE]->empty() )
        {
          parseParam_();
          computeFeatureBuckets_();
          computeShiftBuckets_();
          computeShift_();
          // std::cout << "GeomHashShiftSuperimposer::run()  fertig" << std::endl;
        }
        else
        {
          std::cerr << "GeomHashShiftSuperimposer::run():  Oops, one of the feature maps is empty!\n";
        }
      }

    protected:

      //----------------------------------------------------------------------

      /** @name Methods
       */
      //@{

      void parseParam_()
      {
#define V_parseParam_(bla) V_GeomHashShiftSuperimposer(bla)
        V_parseParam_("@@@ parseParam_()");

        // Initialize feature_bucket_size_ with values from param_.
        std::string fm_bs = "feature_map:bucket_size:";
        for ( Size dimension = 0; dimension < 2; ++dimension)
        {
          std::string fm_bs_dn = fm_bs + DimensionDescriptionType::dimension_name_short[dimension];
          DataValue data_value = getParam().getValue(fm_bs_dn);
          if ( data_value == DataValue::EMPTY )
          {
            throw Exception::ElementNotFound<std::string>
            (__FILE__,__LINE__,__PRETTY_FUNCTION__,fm_bs_dn);
          }
          else
          {
            feature_bucket_size_[dimension] = data_value;
            V_parseParam_(fm_bs_dn<< ": "<<feature_bucket_size_[dimension]);
          }
        }

        // Initialize shift_bucket_size_ with values from param_.
        std::string tm_bs  = "shift_map:bucket_size:";
        for ( Size dimension = 0; dimension < 2; ++dimension)
        {
          String tm_bs_dn
          = tm_bs
            + DimensionDescriptionType::dimension_name_short[dimension];
          DataValue data_value = getParam().getValue(tm_bs_dn);
          if ( data_value == DataValue::EMPTY )
          {
            throw Exception::ElementNotFound<std::string>
            (__FILE__,__LINE__,__PRETTY_FUNCTION__,tm_bs_dn);
          }
          else
          {
            shift_bucket_size_[dimension] = data_value;
            V_parseParam_(tm_bs_dn<< ": "<<shift_bucket_size_[dimension]);
          }
        }

        // Initialize feature_bucket_window_ with values from param_.
        std::string tm_fbw = "feature_map:bucket_window:";
        for ( Size dimension = 0; dimension < 2; ++dimension)
        {
          std::string tm_fbw_dn
          = tm_fbw
            + DimensionDescriptionType::dimension_name_short[dimension];
          DataValue data_value = getParam().getValue(tm_fbw_dn);
          if ( data_value == DataValue::EMPTY )
          {
            throw Exception::ElementNotFound<std::string>
            (__FILE__,__LINE__,__PRETTY_FUNCTION__,tm_fbw_dn);
          }
          else
          {
            feature_bucket_window_[dimension] = data_value;
            V_parseParam_(tm_fbw_dn<< ": "<<feature_bucket_window_[dimension]);
          }
        }

        // Initialize shift_bucket_window_ with values from param_.
        std::string const tm_tbw = "shift_map:bucket_window:";
        for ( Size dimension = 0; dimension < 2; ++dimension)
        {
          std::string tm_tbw_dn
          = tm_tbw
            + DimensionDescriptionType::dimension_name_short[dimension];
          DataValue data_value = getParam().getValue(tm_tbw_dn);
          if ( data_value == DataValue::EMPTY )
          {
            throw Exception::ElementNotFound<std::string>
            (__FILE__,__LINE__,__PRETTY_FUNCTION__,tm_tbw_dn);
          }
          else
          {
            shift_bucket_window_[dimension] = data_value;
            V_parseParam_(tm_tbw_dn<< ": "<<shift_bucket_window_[dimension]);
          }
        }

#undef V_parseParam_

      }


      /**@brief Fill the buckets with the indices of the corresponding features.
       */
      void computeFeatureBuckets_()
      {
#define V_computeFeatureBuckets_(bla) V_GeomHashShiftSuperimposer(bla)
        V_computeFeatureBuckets_("@@@ computeFeatureBuckets_()");

        // Shorthands ...
        PositionType & fbs = feature_bucket_size_;

        for ( Size map_index = 0; map_index < 2; ++map_index )
        {
          // Shorthands ...
          V_computeFeatureBuckets_("\n--- map_index: "<<map_index);
          PointMapType const     & fm     = getFeatureMap(map_index);
          PositionBoundingBoxType  & fmpbb  = feature_map_position_bounding_box_[map_index] ;
          IntensityBoundingBoxType & fmibb  = feature_map_intensity_bounding_box_[map_index];

          fmpbb.clear();
          fmibb.clear();

          // Compute the bounding box for the feature map, with respect to
          // position and intensity.
          for ( typename PointMapType::ConstIterator fm_iter = fm.begin();
                fm_iter != fm.end();
                ++fm_iter
              )
          {
            fmpbb.enlarge(fm_iter->getPosition());
            fmibb.enlarge(fm_iter->getIntensity());
          }
          V_computeFeatureBuckets_("fmpbb: "<<fmpbb<<"fmibb: "<<fmibb);
        }

        // Next we will enlarge each feature_map_position_bounding_box_ such
        // that all buckets will have the same diagonal.  To provide against
        // rounding errors, we allocate one bucket more than needed (in each
        // dimension) and shift the grid by one-half of the difference.
        for ( Size map_index = 0; map_index < 2; ++map_index )
        {
          // Shorthands ...
          V_computeFeatureBuckets_("\n--- map_index: "<<map_index);
          PointMapType          const & fm     = getFeatureMap(map_index);
          PositionBoundingBoxType const & fmpbb  = feature_map_position_bounding_box_[map_index] ;
          PositionBoundingBoxType       & fmpbbe = feature_map_position_bounding_box_enlarged_[map_index] ;
          FeatureBucketMatrixType       & fb     = feature_bucket_[map_index];

          // Compute num_buckets.  Compute extra margin to make bounding box a
          // multiple of feature buckets.
          PositionType const diagonal = fmpbb.diagonal();
          PositionType diagonal_enlarged;
          V_computeFeatureBuckets_("diagonal: " << diagonal);
          int num_buckets[2];
          for ( Size dimension = 0; dimension < 2; ++dimension)
          {
            num_buckets[dimension] = int(1.1 + diagonal[dimension]/fbs[dimension]);
            diagonal_enlarged[dimension] = fbs[dimension] * num_buckets[dimension];
          }
          V_computeFeatureBuckets_("num_buckets: "<<num_buckets[RT]<<' '<<num_buckets[MZ]);
          V_computeFeatureBuckets_("diagonal_enlarged: "<<diagonal_enlarged);

          // The extra margin.
          PositionType extra_feature_bucket_size_(diagonal_enlarged-diagonal);
          extra_feature_bucket_size_ /= 2;
          V_computeFeatureBuckets_("efbs: " << extra_feature_bucket_size_);

          // Compute the enlarged feature map bounding box accordingly.
          fmpbbe.clear();
          fmpbbe.enlarge( fmpbb.min() - extra_feature_bucket_size_ );
          fmpbbe.enlarge( fmpbb.max() + extra_feature_bucket_size_ );
          V_computeFeatureBuckets_("fmpbbe: "<<fmpbbe);

          // Resize feature_bucket_[map_index] accordingly.
          fb.resize(num_buckets[RT],num_buckets[MZ]);
          V_computeFeatureBuckets_("rows: "<<fb.rows()<<"  cols: "<<fb.cols());

          // Now, finally, we store the indices of the features in their
          // corresponding buckets.
          PositionType const & fmpbbe_min = fmpbbe.min();
          for ( Size index= 0; index < fm.size(); ++index )
          {
            PositionType position = fm[index].getPosition() - fmpbbe_min;
            fb ( Size(position[RT]/fbs[RT]), Size(position[MZ]/fbs[MZ]) ).push_back(index);
          }

          // Optionally, write debug output as specified in param.
          String feature_buckets_file_base = getParam().getValue("debug:feature_buckets_file");
          if ( !feature_buckets_file_base.empty() )
          {
            String const feature_buckets_file = feature_buckets_file_base+String(map_index?"_SCENE":"_MODEL");
            std::ofstream dump_file(feature_buckets_file.c_str());
            std::cerr << "### Writing "<<feature_buckets_file<<std::endl;
            dump_file << "# " << feature_buckets_file << " generated " << Date::now() << std::endl;
            dump_file << "# Positions of features in non-empty feature buckets" << std::endl;
            for ( FeatureBucketMatrixType::ConstIterator iter = fb.begin(); iter != fb.end(); ++iter)
            {
              if (iter->empty())
                continue;
              std::pair<Size,Size> row_col = fb.indexPair(iter-fb.begin());
              dump_file << row_col.first << ' ' << row_col.second << " #bucket" << std::endl;
              for ( FeatureBucketType::const_iterator viter = iter->begin(); viter != iter->end(); ++viter)
              {
                dump_file << fm[*viter].getPosition()[RT] <<' '<<fm[*viter].getPosition()[MZ] << std::endl;
              }
              dump_file << std::endl;
            }
            dump_file << "# " << feature_buckets_file << " EOF " << Date::now() << std::endl;
          }
        }

        return;
#undef V_computeFeatureBuckets_

      } // computeFeatureBuckets_


      //----------------------------------------------------------------------

      /**@brief Fill the buckets of shifts.

      Note that computeFeatureBuckets_() must have been called before to make
      this work properly.
      */
      void computeShiftBuckets_()
      {
#define V_computeShiftBuckets_(bla) V_GeomHashShiftSuperimposer(bla)
        V_computeShiftBuckets_("\n");
        V_computeShiftBuckets_("@@@ computeShiftBuckets_()");

        // Shorthands ...
        ShiftQualityMatrixType & tb      = shift_bucket_;
        PositionType                 & tbs     = shift_bucket_size_;
        PositionBoundingBoxType      & tbb     = shift_bounding_box_ ;
        PositionBoundingBoxType      & tbbe    = shift_bounding_box_enlarged_ ;
        Size                   const (&fbw)[2] = feature_bucket_window_;
        ShiftMatrixType        & tm      = shift_matrix_;

        // Compute the bounding box for the shift map
        {
          tbb.clear();
          tbb.enlarge ( feature_map_position_bounding_box_[MZ].min() - feature_map_position_bounding_box_[RT].min() );
          tbb.enlarge ( feature_map_position_bounding_box_[MZ].min() - feature_map_position_bounding_box_[RT].max() );
          tbb.enlarge ( feature_map_position_bounding_box_[MZ].max() - feature_map_position_bounding_box_[RT].min() );
          tbb.enlarge ( feature_map_position_bounding_box_[MZ].max() - feature_map_position_bounding_box_[RT].max() );
        }
        V_computeShiftBuckets_("tbb: "<<tbb);

        // Next we will enlarge each bucket_size_ such that all buckets will
        // have the same diagonal.  To provide against rounding errors, we
        // allocate one bucket more than needed (in each dimension) and shift
        // the grid by one-half.

        PositionType half_of_shift_bucket_size_(tbs);
        half_of_shift_bucket_size_ /= 2;
        V_computeShiftBuckets_("hotbs: " << half_of_shift_bucket_size_);

        // Adjust the enlarged shift map bounding box accordingly.
        {
          tbbe.clear();
          tbbe.enlarge( tbb.min() - half_of_shift_bucket_size_ );
          tbbe.enlarge( tbb.max() + half_of_shift_bucket_size_ );
        }
        V_computeShiftBuckets_("tbbe: "<<tbbe);

        // Compute shift_bucket_size_ and num_buckets.
        PositionType diagonal = tbbe.diagonal();
        V_computeShiftBuckets_("diagonal: " << diagonal);
        int num_buckets[2];
        for ( Size dimension = 0; dimension < 2; ++dimension)
        {
          num_buckets[dimension] = int(diagonal[dimension]/tbs[dimension]);
          tbs[dimension] = diagonal[dimension] / num_buckets[dimension];
        }
        V_computeShiftBuckets_("tbs: "<<tbs);

        // Resize shift_bucket_ accordingly.
        tb.resize(num_buckets[RT]+1,num_buckets[MZ]+1);
        V_computeShiftBuckets_("rows: "<<tb.rows()<<"  cols: "<<tb.cols());

        // Clear the shift buckets.
        std::fill(tb.begin(),tb.end(),QualityType(0));


        // Resize shift_matrix_ according to feature_bucket_[MZ]
        tm.resize(feature_bucket_[MZ].sizePair());

        // Now we store the shifts for all relevant feature pairs in their
        // corresponding buckets.  Each shift is distributed among its
        // four neighboring "buckets", with weights according to the distances
        // from these corner points.  Note that the outer two loops (over i and
        // j) enumerate the "image" (feature_bucket_[MZ]), then we search for
        // "pre-images" (feature_bucket_[0}) in the two inner loops (over k and
        // l).  (And of course, finally, we enumerate all feature pairs.)  This
        // way we can associate the shifts vectors to buckets of the
        // image, and when we will later apply it, we will not change the
        // pre-image, which might be a consensus or so.
#define V_computeShiftBuckets_enumeration(bla) V_computeShiftBuckets_(bla)

        // progress dots
        DataValue const & param_progress_dots = this->param_.getValue("debug:progress_dots");
        int progress_dots
        = param_progress_dots.isEmpty() ? 0 : int(param_progress_dots);

        PositionType const & tbbe_min = tbbe.min();

        // Compute the index shift of corresponding feature buckets of model and scene.
        PositionType const fmpbbe_min_offset =
          feature_map_position_bounding_box_enlarged_[SCENE].min() -
          feature_map_position_bounding_box_enlarged_[MODEL].min();
        int const feature_buckets_index_offset_RT = int ( fmpbbe_min_offset[RT] / feature_bucket_size_[RT] );
        int const feature_buckets_index_offset_MZ = int ( fmpbbe_min_offset[MZ] / feature_bucket_size_[MZ] );

        // iterate over buckets of scene
        for ( Size scene_bucket_index_RT = 0;
              scene_bucket_index_RT < feature_bucket_[SCENE].rows();
              ++scene_bucket_index_RT
            )
        {
          for ( Size scene_bucket_index_MZ = 0;
                scene_bucket_index_MZ < feature_bucket_[SCENE].cols();
                ++scene_bucket_index_MZ
              )
          {

            // compute the corresponding bucket in the model
            int const model_bucket_index_center_RT = scene_bucket_index_RT + feature_buckets_index_offset_RT;
            int const model_bucket_index_center_MZ = scene_bucket_index_MZ + feature_buckets_index_offset_MZ;

            // iterate over buckets of model
            for ( int model_bucket_index_RT
                  =  std::max<int>( model_bucket_index_center_RT - fbw[RT], 0 );
                  model_bucket_index_RT
                  <= std::min<int>( model_bucket_index_center_RT + fbw[RT], feature_bucket_[MODEL].rows()-1 );
                  ++model_bucket_index_RT
                )
            {
              for ( int model_bucket_index_MZ
                    =  std::max<int>( model_bucket_index_center_MZ - fbw[MZ], 0 );
                    model_bucket_index_MZ
                    <= std::min<int>( model_bucket_index_center_MZ + fbw[MZ], feature_bucket_[MODEL].cols()-1 );
                    ++model_bucket_index_MZ
                  )
              {
                // iterate over pairs of features for this pair of buckets
                int number_of_considered_feature_pairs_for_this_pair_of_buckets = 0;
                FeatureBucketType const & model_feature_bucket
                = feature_bucket_[MODEL]
                  ( model_bucket_index_RT, model_bucket_index_MZ );
                for ( FeatureBucketType::const_iterator model_iter = model_feature_bucket.begin();
                      model_iter != model_feature_bucket.end();
                      ++model_iter
                    )
                {
                  FeatureBucketType const & scene_feature_bucket
                  = feature_bucket_[SCENE]( scene_bucket_index_RT, scene_bucket_index_MZ );
                  for ( FeatureBucketType::const_iterator scene_iter = scene_feature_bucket.begin();
                        scene_iter != scene_feature_bucket.end();
                        ++scene_iter
                      )
                  {
                    // Compute the shift corresponding to a pair of features.
                    ShiftType shift = shift_( getFeatureMap(0)[*model_iter],
                                              getFeatureMap(1)[*scene_iter] );
                    //                     V_computeShiftBuckets_enumeration("shift: "<< shift.getPosition());
                    //                     V_computeShiftBuckets_enumeration("shift: "<< shift.getQuality());

                    PositionType tpwm = shift.getPosition();
                    tpwm -= tbbe_min;
                    // V_computeShiftBuckets_enumeration("trans.pos wrt tbbe_min: "<< shift.getPosition());

                    QualityType  const & tq = shift.getQuality();

                    // Compute the bucket index (the lowest of the four) for
                    // this shift.  Also compute the fractional part of
                    // the position within the bucket.
                    Size bucket_index[2];
                    PositionType bucket_fraction;
                    for ( Size dimension = 0; dimension < 2; ++dimension )
                    {
                      bucket_fraction[dimension] = tpwm[dimension] / tbs[dimension];  // floating point division
                      bucket_index[dimension]    = (Size) bucket_fraction[dimension]; // round down (yes we are >= 0)
                      bucket_fraction[dimension] -= bucket_index[dimension];          // fractional part
                    }
                    PositionType bucket_fraction_complement(1,1);
                    bucket_fraction_complement -= bucket_fraction;

                    // Distribute the quality of the shift among the four neighboring buckets.
                    QualityType factor;

                    factor = bucket_fraction_complement[RT] * bucket_fraction_complement[MZ];
                    tb( bucket_index[RT], bucket_index[MZ] ) += tq * factor;

                    factor = bucket_fraction_complement[RT] * bucket_fraction[MZ];
                    tb( bucket_index[RT], bucket_index[MZ] + 1 ) += tq * factor;

                    factor = bucket_fraction[RT] * bucket_fraction_complement[MZ];
                    tb( bucket_index[RT] + 1, bucket_index[MZ] ) += tq * factor;

                    factor = bucket_fraction[RT] * bucket_fraction[MZ];
                    tb( bucket_index[RT] + 1, bucket_index[MZ] + 1 ) += tq * factor;

                    ++number_of_considered_feature_pairs_for_this_pair_of_buckets;

                    if ( progress_dots &&
                         ! (number_of_considered_feature_pairs_for_this_pair_of_buckets % progress_dots)
                       )
                    {
											std::cout << 'H' << std::flush;
                    }

                  } // for scene_iter
                } // for model_iter

#if 0 // debug output
                if ( number_of_considered_feature_pairs_for_this_pair_of_buckets )
                {
                  std::cout <<
                  "s_b_i_RT, _MZ, m_b_i_c_RT, _MZ, m_b_i_RT, _MZ, number_pairs: " <<
                  scene_bucket_index_RT<<' '<<scene_bucket_index_MZ<<' '<<
                  model_bucket_index_center_RT<<' '<<model_bucket_index_center_MZ<<' '<<
                  model_bucket_index_RT<<' '<<model_bucket_index_MZ<<' '<<
                  number_of_considered_feature_pairs_for_this_pair_of_buckets
                  ;
                }
#endif

              } // for model_bucket_index_MZ
            } // for model_bucket_index_RT
          } // for scene_bucket_index_MZ
        } // for scene_bucket_index_RT

#undef V_computeShiftBuckets_enumeration

        // Optionally, write debug output as specified in param.
        DataValue data_value_dump_shift_buckets = getParam().getValue("debug:dump_shift_buckets");
        if ( data_value_dump_shift_buckets != DataValue::EMPTY )
        {
          std::string   dump_filename = data_value_dump_shift_buckets;
          std::ofstream dump_file(dump_filename.c_str());
          V_computeShiftBuckets_("### Writing "<<dump_filename);
          dump_file << "# " << dump_filename << " generated " << Date::now() << std::endl;
          dump_file << "# Shift buckets: xcoord ycoord quality xindex yindex" << std::endl;

          for ( typename ShiftQualityMatrixType::ConstIterator iter = tb.begin(); iter != tb.end(); ++iter)
          {
            std::pair<Size,Size> row_col = tb.indexPair(iter-tb.begin());
            if ( *iter )
            {
              dump_file << tbbe_min[RT] + tbs[RT] * row_col.first << ' '
              << tbbe_min[MZ] + tbs[MZ] * row_col.second << ' '
              << *iter << ' '
              << row_col.first << ' '
              << row_col.second
              << " #tb" << std::endl ;
            }
          }
          dump_file << "# " << dump_filename << " EOF " << Date::now() << std::endl;
        }

#undef V_computeShiftBuckets_

      } // computeShiftBuckets_



      //----------------------------------------------------------------------

      /**@brief Compute the shift.

      Note that shift_buckets_ must have been calculated before.
      */
      void computeShift_()
      {
#define V_computeShift_(bla) V_GeomHashShiftSuperimposer(bla)
        V_computeShift_("@@@ computeShift_()");

        ShiftType shift;

        // Shorthands ...
        ShiftQualityMatrixType const & tb = shift_bucket_;
        PositionType const & tbs = shift_bucket_size_;
        Size const (&tbw)[2] = shift_bucket_window_;

        // Find the transformation bucket with highest impact (quality).
        Size tb_max_element_index = std::max_element(tb.begin(),tb.end()) - tb.begin();
        Size tb_max_indices[2];
        tb_max_indices[RT] = tb.rowIndex(tb_max_element_index);
        tb_max_indices[MZ] = tb.colIndex(tb_max_element_index);
        V_computeShift_("tb_max: "<<tb_max_indices[RT]<<' '<<tb_max_indices[MZ]<<" quality="<<tb(tb_max_indices[RT],tb_max_indices[MZ]));

        // Compute a weighted average of the shifts nearby the tb_max_element.
        //ShiftType result; // initially zero

        PositionType const& tbbe_min = shift_bounding_box_enlarged_.min();
        int tb_run_indices[2];
        for ( tb_run_indices[RT]  = std::max ( int (tb_max_indices[RT] - tbw[RT]), 0 );
              tb_run_indices[RT] <= std::min ( int (tb_max_indices[RT] + tbw[RT]), int (tb.rows()) - 1 );
              ++tb_run_indices[RT]
            )
        {
          for ( tb_run_indices[MZ]  = std::max ( int (tb_max_indices[MZ] - tbw[MZ]), 0 );
                tb_run_indices[MZ] <= std::min ( int (tb_max_indices[MZ] + tbw[MZ]), int (tb.cols()) - 1 );
                ++tb_run_indices[MZ]
              )
          {
            PositionType contribution_position(tbs);
            for ( Size dimension = 0; dimension < 2; ++dimension)
            {
              contribution_position[dimension] *= tb_run_indices[dimension];
            }
            contribution_position += tbbe_min;
            QualityType contribution_quality = tb( tb_run_indices[RT], tb_run_indices[MZ] );
            shift.getQuality() += contribution_quality;
            contribution_position *= contribution_quality;
            shift.getPosition() += contribution_position;
          }
        }
        if ( shift.getQuality() != 0 )
        {
          // ???? we need to shift the opposite way!  Found that out by try-and-error.  Clemens
          shift.getPosition() /= -shift.getQuality() ;
        }
        else
        {
          // result.getPosition() is irrelevant anyway
        }

        // Assign the result.
        for ( int dim = 0; dim < 2; ++dim )
        {
          if ( !final_transformation_[dim] )
          {
            final_transformation_[dim] = new FinalShiftType;
          }
          FinalShiftType *f = dynamic_cast<FinalShiftType*>( final_transformation_[dim] );
          // set slope and intercept
          f->setParam( 1.0, shift.getPosition()[dim] );
          V_computeShift_("computeShift_() hat geklappt: " << shift.getPosition());
        }



#undef V_computeShift_

      } // computeShift_

      //----------------------------------------------------------------------

      /**@brief Compute the shift and similarity for a pair of features;
         larger quality values are better.

         The returned value should express our confidence that one feature might
         possibly be matched to the other.

         Currently this will just calculate the ratio of intensities, either
         "left/right" or "right/left", such that a value between 0 and 1 is
         returned.

         \todo Take the quality of the features themselves into account, i.e.,
         how good they fit to their model.
      */
      ShiftType shift_( PointType const & left, PointType const & right ) const
      {

        // @todo Take the quality of the features themselves into account, i.e. how good they fit to their model.

        ShiftType shift;
        shift.setPosition(right.getPosition() - left.getPosition());
        if ( right.getIntensity() == 0 )
          shift.setQuality(0);
        QualityType result = left.getIntensity() / right.getIntensity();
        shift.setQuality( result <= 1. ? result : 1. / result );
        return shift;
      }

      //@} // Methods


      /** @name Data members
       */
      //@{

      /// Holds the bounding box of all input features.
      PositionBoundingBoxType  feature_map_position_bounding_box_[2];

      /// Holds the enlarged bounding box for all input features.  It is larger
      /// by about half of a bucket in all directions.
      PositionBoundingBoxType  feature_map_position_bounding_box_enlarged_[2];

      /// Holds a bounding box for the input feature intensities.
      IntensityBoundingBoxType feature_map_intensity_bounding_box_[2];

      /// Feature indices are stored in theses buckets.
      FeatureBucketMatrixType feature_bucket_[2];

      /// Diagonal size of each bucket.
      PositionType feature_bucket_size_;

      /// Shifts are stored (summed up) in these buckets.
      ShiftQualityMatrixType shift_bucket_;

      /// Holds a bounding box for all possible shift vectors.
      PositionBoundingBoxType shift_bounding_box_;

      /// Holds an enlarged bounding box for all shift vectors.  It is
      /// larger by about half of a bucket in all directions.
      PositionBoundingBoxType shift_bounding_box_enlarged_;

      /// Diagonal size of each bucket in shift_bucket_.
      PositionType shift_bucket_size_;

      /// Number of surrounding buckets of feature indices to be considered when
      /// computing shifts.
      Size feature_bucket_window_[2];

      /// Number of surrounding buckets of shift indices to be considered when
      /// computing shifts.
      Size shift_bucket_window_[2];

      /// Matrix of shifts associated with buckets of feature_map_[SCENE].
      ShiftMatrixType shift_matrix_;
      //@}

  }
  ; // GeomHashShiftSuperimposer

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_GeomHashShiftSuperimposer_H
