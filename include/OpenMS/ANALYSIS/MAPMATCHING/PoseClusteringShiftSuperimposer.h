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
// $Maintainer: Clemens Groepl, Eva Lange $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGSHIFTSUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGSHIFTSUPERIMPOSER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DLinearMapping.h>


#include <fstream>

#define V_PoseClusteringShiftSuperimposer(bla) // std::cout << bla << std::endl;

namespace OpenMS
{

  /**
  @brief Superimposer that uses a voting scheme to find a good translation.

     It works on two element maps (DFeatureMap is the default map type, 
     but you can also use a pointer map like DPeakConstReferenceArray) and 
     computes a translation, that maps the elements of one map (scene map) 
     as near as possible to the elements in the other map (model map).
     A element can be a DPeak, a DFeature, a ConsensusPeak or ConsensusFeature 
     (wheras DFeature is the default element type).
     
     This superimposer hashs all possible shifts and defines the 
     translation with the most votes as the best one.        
  */
  template < typename MapT = DFeatureMap<2> >
  class PoseClusteringShiftSuperimposer
        : public BaseSuperimposer< MapT >
  {
  public:
    /// Defines the coordinates of peaks / elements.
    typedef DimensionDescription<LCMS_Tag> DimensionDescriptionType;
    enum DimensionId
    {
      RT = DimensionDescriptionType::RT,
      MZ = DimensionDescriptionType::MZ
    };

    /** Symbolic names for indices of element maps etc.
        This should make things more understandable and maintainable.
         */
    enum Maps
    {
      MODEL = 0,
      SCENE = 1
    };

    typedef BaseSuperimposer< MapT > Base;

  protected:
    /**
    @brief Intensity bounding box

    We need this to make the intensity bounding box use the intensity type
    instead of the coordinate type.
    */
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

    typedef typename Base::TraitsType TraitsType;
    typedef typename Base::QualityType QualityType;
    typedef typename Base::PositionType PositionType;
    typedef typename Base::IntensityType IntensityType;
    typedef typename Base::PointType PointType;
    typedef typename Base::PointMapType PointMapType;
    typedef typename PositionType::CoordinateType CoordinateType;
    typedef DBoundingBox<2,TraitsType>  PositionBoundingBoxType;
    typedef DBoundingBox<1,IntensityBoundingBoxTraits> IntensityBoundingBoxType;
    typedef std::vector <Size> ElementBucketType;
    typedef Matrix < ElementBucketType > ElementBucketMatrixType;
    typedef Shift ShiftType;
    typedef Matrix < typename ShiftType::QualityType > ShiftQualityMatrixType;
    typedef Matrix < ShiftType > ShiftMatrixType;
    typedef DLinearMapping< 1, TraitsType > FinalShiftType;

    using Base::setParameters;
    using Base::getParameters;
    using Base::subsections_;
    using Base::defaultsToParam_;
    using Base::param_;
    using Base::defaults_;
    using Base::setElementMap;
    using Base::getElementMap;
    using Base::final_transformation_;

    /// Constructor
    PoseClusteringShiftSuperimposer()
        : Base()
    {
			setName(getProductName());
			
      defaults_.setValue("feature_map:bucket_size:RT",150);
      defaults_.setValue("feature_map:bucket_size:MZ",4);
      defaults_.setValue("transformation_space:shift_bucket_size:RT",5);
      defaults_.setValue("transformation_space:shift_bucket_size:MZ",0.1);
      defaults_.setValue("feature_map:bucket_window:RT",2);
      defaults_.setValue("feature_map:bucket_window:MZ",1);
      defaults_.setValue("transformation_space:bucket_window_shift:RT",2);
      defaults_.setValue("transformation_space:bucket_window_shift:MZ",1);
			subsections_.push_back("debug");
			
      defaultsToParam_();
    }

    /// Copy constructor
    PoseClusteringShiftSuperimposer(const PoseClusteringShiftSuperimposer& source)
        : Base(source),
        element_bucket_(source.element_bucket_),
        shift_bucket_(source.shift_bucket_),
        shift_bounding_box_(source.shift_bounding_box_),
        shift_bounding_box_enlarged_(source.shift_bounding_box_enlarged_)
    {
			updateMembers_();
    }

    ///  Assignment operator
    PoseClusteringShiftSuperimposer& operator = (const PoseClusteringShiftSuperimposer& source)
    {
      if (&source==this) return *this;

      Base::operator=(source);
      	
      element_bucket_[0] = source.element_bucket_[0];
      element_bucket_[1] = source.element_bucket_[1];
      shift_bucket_ = source.shift_bucket_;
      shift_bounding_box_ = source.shift_bounding_box_;
      shift_bounding_box_enlarged_ = source.shift_bounding_box_enlarged_;
 			
 			updateMembers_();
 			
      return *this;
    }

    /// Destructor
    virtual ~PoseClusteringShiftSuperimposer()
    {
      V_PoseClusteringShiftSuperimposer("~PoseClusteringShiftSuperimposer");
    }

    /// Estimates the transformation for each grid cell
    virtual void run()
    {
      // clear the member
      element_bucket_[RT].clear();
      element_bucket_[MZ].clear();
      shift_bucket_.clear();

      if ( !this->element_map_[MODEL]->empty() && !this->element_map_[SCENE]->empty() )
      {
        computeElementBuckets_();
        computeShiftBuckets_();
        computeShift_();
      }
      else
      {
        std::cerr << "PoseClusteringShiftSuperimposer::run():  Oops, one of the element maps is empty!\n";
      }
    }

    /// Returns an instance of this class
    static BaseSuperimposer<PointMapType>* create()
    {
      return new PoseClusteringShiftSuperimposer();
    }

    /// Returns the name of this module
    static const String getProductName()
    {
      return "poseclustering_shift";
    }

    /// Set size of shift buckets (in dimension dim)
    void setShiftBucketSize(UnsignedInt dim, double shift_bucket_size)
    {
      shift_bucket_size_[dim] = shift_bucket_size;
      param_.setValue( String("transformation_space:shift_bucket_size:") + DimensionDescriptionType::dimension_name_short[dim], (float)shift_bucket_size);
    }

    /// Get size of shift buckets (in dimension dim)
    double getShiftBucketSize(UnsignedInt dim) const
    {
      return shift_bucket_size_[dim];
    }

    /// Set number of neighbouring element buckets to be considered for the calculation of the final transformation (in dimension dim)
    void setElementBucketWindow(UnsignedInt dim, Size element_bucket_window)
    {
      element_bucket_window_[dim] = element_bucket_window;
      param_.setValue(String("feature_map:bucket_window:") + DimensionDescriptionType::dimension_name_short[dim], (int)element_bucket_window);
    }

    /// Get number of neighbouring shift buckets to be considered for the calculation of the final transformation (in dimension dim)
    Size getElementBucketWindow(UnsignedInt dim) const
    {
      return element_bucket_window_[dim];
    }

    /// Set number of neighbouring shift buckets to be considered for the calculation of the final transformation (in dimension dim)
    void setShiftBucketWindow(UnsignedInt dim, Size shift_bucket_window)
    {
      shift_bucket_window_[dim] = shift_bucket_window;
      param_.setValue(String("transformation_space:bucket_window_shift:") + DimensionDescriptionType::dimension_name_short[dim], (int)shift_bucket_window);
    }

    /// Get number of neighbouring shift buckets to be considered for the calculation of the final transformation (in dimension dim)
    Size getShiftBucketWindow(UnsignedInt dim) const
    {
      return shift_bucket_window_[dim];
    }
  protected:
    virtual void updateMembers_()
    {
      shift_bucket_size_[0] = (CoordinateType)param_.getValue("transformation_space:shift_bucket_size:RT");
      shift_bucket_size_[1] = (CoordinateType)param_.getValue("transformation_space:shift_bucket_size:MZ");
      element_bucket_window_[0] = (Size)param_.getValue("feature_map:bucket_window:RT");
      element_bucket_window_[1] = (Size)param_.getValue("feature_map:bucket_window:MZ");
      shift_bucket_window_[0] = (Size)param_.getValue("transformation_space:bucket_window_shift:RT");
      shift_bucket_window_[1] = (Size)param_.getValue("transformation_space:bucket_window_shift:MZ");
      element_bucket_size_[0] = (CoordinateType)param_.getValue("feature_map:bucket_size:RT");
      element_bucket_size_[1] = (CoordinateType)param_.getValue("feature_map:bucket_size:MZ");
    }

    /// Fill the buckets with the indices of the corresponding elements.
    void computeElementBuckets_()
    {
#define V_computeElementBuckets_(bla) V_PoseClusteringShiftSuperimposer(bla)
      V_computeElementBuckets_("@@@ computeElementBuckets_()");

      // Shorthands ...
      PositionType & fbs = element_bucket_size_;

      for ( Size map_index = 0; map_index < 2; ++map_index )
      {
        // Shorthands ...
        V_computeElementBuckets_("\n--- map_index: "<<map_index);
        PointMapType const     & fm     = getElementMap(map_index);
        PositionBoundingBoxType  & fmpbb  = element_map_position_bounding_box_[map_index] ;
        IntensityBoundingBoxType & fmibb  = element_map_intensity_bounding_box_[map_index];

        fmpbb.clear();
        fmibb.clear();

        // Compute the bounding box for the element map, with respect to
        // position and intensity.
        for ( typename PointMapType::ConstIterator fm_iter = fm.begin();
              fm_iter != fm.end();
              ++fm_iter
            )
        {
          fmpbb.enlarge(fm_iter->getPosition());
          fmibb.enlarge(fm_iter->getIntensity());
        }
        V_computeElementBuckets_("fmpbb: "<<fmpbb<<"fmibb: "<<fmibb);
      }

      // Next we will enlarge each element_map_position_bounding_box_ such
      // that all buckets will have the same diagonal.  To provide against
      // rounding errors, we allocate one bucket more than needed (in each
      // dimension) and shift the grid by one-half of the difference.
      for ( Size map_index = 0; map_index < 2; ++map_index )
      {
        // Shorthands ...
        V_computeElementBuckets_("\n--- map_index: "<<map_index);
        PointMapType          const & fm     = getElementMap(map_index);
        PositionBoundingBoxType const & fmpbb  = element_map_position_bounding_box_[map_index] ;
        PositionBoundingBoxType       & fmpbbe = element_map_position_bounding_box_enlarged_[map_index] ;
        ElementBucketMatrixType       & fb     = element_bucket_[map_index];

        // Compute num_buckets.  Compute extra margin to make bounding box a
        // multiple of element buckets.
        PositionType const diagonal = fmpbb.diagonal();
        PositionType diagonal_enlarged;
        V_computeElementBuckets_("diagonal: " << diagonal);
        int num_buckets[2];
        for ( Size dimension = 0; dimension < 2; ++dimension)
        {
          num_buckets[dimension] = int(1.1 + diagonal[dimension]/fbs[dimension]);
          diagonal_enlarged[dimension] = fbs[dimension] * num_buckets[dimension];
        }
        V_computeElementBuckets_("num_buckets: "<<num_buckets[RT]<<' '<<num_buckets[MZ]);
        V_computeElementBuckets_("diagonal_enlarged: "<<diagonal_enlarged);

        // The extra margin.
        PositionType extra_element_bucket_size_(diagonal_enlarged-diagonal);
        extra_element_bucket_size_ /= 2;
        V_computeElementBuckets_("efbs: " << extra_element_bucket_size_);

        // Compute the enlarged element map bounding box accordingly.
        fmpbbe.clear();
        fmpbbe.enlarge( fmpbb.min() - extra_element_bucket_size_ );
        fmpbbe.enlarge( fmpbb.max() + extra_element_bucket_size_ );
        V_computeElementBuckets_("fmpbbe: "<<fmpbbe);

        // Resize element_bucket_[map_index] accordingly.
        fb.resize(num_buckets[RT],num_buckets[MZ]);
        V_computeElementBuckets_("rows: "<<fb.rows()<<"  cols: "<<fb.cols());

        // Now, finally, we store the indices of the elements in their
        // corresponding buckets.
        PositionType const & fmpbbe_min = fmpbbe.min();
        for ( Size index= 0; index < fm.size(); ++index )
        {
          PositionType position = fm[index].getPosition() - fmpbbe_min;
          fb ( Size(position[RT]/fbs[RT]), Size(position[MZ]/fbs[MZ]) ).push_back(index);
        }

        // Optionally, write debug output as specified in param.
        String element_buckets_file_base = getParameters().getValue("debug:feature_buckets_file");
        if ( !element_buckets_file_base.empty() )
        {
          String const element_buckets_file = element_buckets_file_base+String(map_index?"_SCENE":"_MODEL");
          std::ofstream dump_file(element_buckets_file.c_str());
          std::cerr << "### Writing "<<element_buckets_file<<std::endl;
          dump_file << "# " << element_buckets_file << " generated " << Date::now() << std::endl;
          dump_file << "# Positions of elements in non-empty element buckets" << std::endl;
          for ( ElementBucketMatrixType::ConstIterator iter = fb.begin(); iter != fb.end(); ++iter)
          {
            if (iter->empty())
              continue;
            std::pair<Size,Size> row_col = fb.indexPair(iter-fb.begin());
            dump_file << row_col.first << ' ' << row_col.second << " #bucket" << std::endl;
            for ( ElementBucketType::const_iterator viter = iter->begin(); viter != iter->end(); ++viter)
            {
              dump_file << fm[*viter].getPosition()[RT] <<' '<<fm[*viter].getPosition()[MZ] << std::endl;
            }
            dump_file << std::endl;
          }
          dump_file << "# " << element_buckets_file << " EOF " << Date::now() << std::endl;
        }
      }

      return;
#undef V_computeElementBuckets_

    } // computeElementBuckets_


    /**@brief Fill the buckets of shifts.

    Note that computeElementBuckets_() must have been called before to make
    this work properly.
    */
    void computeShiftBuckets_()
    {
#define V_computeShiftBuckets_(bla) V_PoseClusteringShiftSuperimposer(bla)
      V_computeShiftBuckets_("\n");
      V_computeShiftBuckets_("@@@ computeShiftBuckets_()");

      // Shorthands ...
      ShiftQualityMatrixType & tb      = shift_bucket_;
      PositionType                 & tbs     = shift_bucket_size_;
      PositionBoundingBoxType      & tbb     = shift_bounding_box_ ;
      PositionBoundingBoxType      & tbbe    = shift_bounding_box_enlarged_ ;
      Size                   const (&fbw)[2] = element_bucket_window_;
      //         ShiftMatrixType        & tm      = shift_matrix_;

      // Compute the bounding box for the shift map
      {
        tbb.clear();
        tbb.enlarge ( element_map_position_bounding_box_[SCENE].min() - element_map_position_bounding_box_[MODEL].min() );
        tbb.enlarge ( element_map_position_bounding_box_[SCENE].min() - element_map_position_bounding_box_[MODEL].max() );
        tbb.enlarge ( element_map_position_bounding_box_[SCENE].max() - element_map_position_bounding_box_[MODEL].min() );
        tbb.enlarge ( element_map_position_bounding_box_[SCENE].max() - element_map_position_bounding_box_[MODEL].max() );
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


      // Resize shift_matrix_ according to element_bucket_[MZ]
      //         tm.resize(element_bucket_[MZ].sizePair());

      // Now we store the shifts for all relevant element pairs in their
      // corresponding buckets.  Each shift is distributed among its
      // four neighboring "buckets", with weights according to the distances
      // from these corner points.  Note that the outer two loops (over i and
      // j) enumerate the "image" (element_bucket_[MZ]), then we search for
      // "pre-images" (element_bucket_[0}) in the two inner loops (over k and
      // l).  (And of course, finally, we enumerate all element pairs.)  This
      // way we can associate the shifts vectors to buckets of the
      // image, and when we will later apply it, we will not change the
      // pre-image, which might be a consensus or so.
#define V_computeShiftBuckets_enumeration(bla) V_computeShiftBuckets_(bla)

      // progress dots
      DataValue const & param_progress_dots = this->getParameters().getValue("debug:progress_dots");
      int progress_dots
      = param_progress_dots.isEmpty() ? 0 : int(param_progress_dots);

      PositionType const & tbbe_min = tbbe.min();

      // Compute the index shift of corresponding element buckets of model and scene.
      PositionType const fmpbbe_min_offset =
        element_map_position_bounding_box_enlarged_[SCENE].min() -
        element_map_position_bounding_box_enlarged_[MODEL].min();
      int const element_buckets_index_offset_RT = int ( fmpbbe_min_offset[RT] / element_bucket_size_[RT] );
      int const element_buckets_index_offset_MZ = int ( fmpbbe_min_offset[MZ] / element_bucket_size_[MZ] );

      // iterate over buckets of scene
      for ( Size scene_bucket_index_RT = 0;
            scene_bucket_index_RT < element_bucket_[SCENE].rows();
            ++scene_bucket_index_RT
          )
      {
        for ( Size scene_bucket_index_MZ = 0;
              scene_bucket_index_MZ < element_bucket_[SCENE].cols();
              ++scene_bucket_index_MZ
            )
        {

          // compute the corresponding bucket in the model
          int const model_bucket_index_center_RT = scene_bucket_index_RT + element_buckets_index_offset_RT;
          int const model_bucket_index_center_MZ = scene_bucket_index_MZ + element_buckets_index_offset_MZ;

          // iterate over buckets of model
          for ( int model_bucket_index_RT
                =  std::max<int>( model_bucket_index_center_RT - fbw[RT], 0 );
                model_bucket_index_RT
                <= std::min<int>( model_bucket_index_center_RT + fbw[RT], element_bucket_[MODEL].rows()-1 );
                ++model_bucket_index_RT
              )
          {
            for ( int model_bucket_index_MZ
                  =  std::max<int>( model_bucket_index_center_MZ - fbw[MZ], 0 );
                  model_bucket_index_MZ
                  <= std::min<int>( model_bucket_index_center_MZ + fbw[MZ], element_bucket_[MODEL].cols()-1 );
                  ++model_bucket_index_MZ
                )
            {
              // iterate over pairs of elements for this pair of buckets
              int number_of_considered_element_pairs_for_this_pair_of_buckets = 0;
              ElementBucketType const & model_element_bucket
              = element_bucket_[MODEL]
                ( model_bucket_index_RT, model_bucket_index_MZ );
              for ( ElementBucketType::const_iterator model_iter = model_element_bucket.begin();
                    model_iter != model_element_bucket.end();
                    ++model_iter
                  )
              {
                ElementBucketType const & scene_element_bucket
                = element_bucket_[SCENE]( scene_bucket_index_RT, scene_bucket_index_MZ );
                for ( ElementBucketType::const_iterator scene_iter = scene_element_bucket.begin();
                      scene_iter != scene_element_bucket.end();
                      ++scene_iter
                    )
                {
                  // Compute the shift corresponding to a pair of elements.
                  ShiftType shift = shift_( getElementMap(0)[*model_iter],
                                            getElementMap(1)[*scene_iter] );
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

                  ++number_of_considered_element_pairs_for_this_pair_of_buckets;

                  if ( progress_dots &&
                       ! (number_of_considered_element_pairs_for_this_pair_of_buckets % progress_dots)
                     )
                  {
                    std::cout << 'H' << std::flush;
                  }

                } // for scene_iter
              } // for model_iter

#if 0 // debug output
              if ( number_of_considered_element_pairs_for_this_pair_of_buckets )
              {
                std::cout <<
                "s_b_i_RT, _MZ, m_b_i_c_RT, _MZ, m_b_i_RT, _MZ, number_pairs: " <<
                scene_bucket_index_RT<<' '<<scene_bucket_index_MZ<<' '<<
                model_bucket_index_center_RT<<' '<<model_bucket_index_center_MZ<<' '<<
                model_bucket_index_RT<<' '<<model_bucket_index_MZ<<' '<<
                number_of_considered_element_pairs_for_this_pair_of_buckets
                ;
              }
#endif

            } // for model_bucket_index_MZ
          } // for model_bucket_index_RT
        } // for scene_bucket_index_MZ
      } // for scene_bucket_index_RT

#undef V_computeShiftBuckets_enumeration

      // Optionally, write debug output as specified in param.
      DataValue data_value_dump_shift_buckets = getParameters().getValue("debug:dump_shift_buckets");
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



    /**@brief Compute the shift.

    Note that shift_buckets_ must have been calculated before.
    */
    void computeShift_()
    {
#define V_computeShift_(bla) V_PoseClusteringShiftSuperimposer(bla)
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
        // set slope and intercept
        final_transformation_[dim].setParam( 1.0, shift.getPosition()[dim] );
        V_computeShift_("computeShift_() hat geklappt: " << shift.getPosition()[dim]);
      }

#undef V_computeShift_

    } // computeShift_


    /**@brief Compute the shift and similarity for a pair of elements;
       larger quality values are better.

       The returned value should express our confidence that one element might
       possibly be matched to the other.

       Currently this will just calculate the ratio of intensities, either
       "left/right" or "right/left", such that a value between 0 and 1 is
       returned.

       @todo Take the quality of the elements themselves into account, i.e., how good they fit to their model.
    */
    ShiftType shift_( PointType const & left, PointType const & right ) const
    {
      ShiftType shift;
      shift.setPosition(right.getPosition() - left.getPosition());
      if ( right.getIntensity() == 0 )
        shift.setQuality(0);
      QualityType result = left.getIntensity() / right.getIntensity();
      shift.setQuality( result <= 1. ? result : 1. / result );
      return shift;
    }

    /// Holds the bounding box of all input elements.
    PositionBoundingBoxType  element_map_position_bounding_box_[2];

    /// Holds the enlarged bounding box for all input elements.  It is larger
    /// by about half of a bucket in all directions.
    PositionBoundingBoxType  element_map_position_bounding_box_enlarged_[2];

    /// Holds a bounding box for the input element intensities.
    IntensityBoundingBoxType element_map_intensity_bounding_box_[2];

    /// Element indices are stored in theses buckets.
    ElementBucketMatrixType element_bucket_[2];

    /// Diagonal size of each bucket.
    PositionType element_bucket_size_;

    /// Shifts are stored (summed up) in these buckets.
    ShiftQualityMatrixType shift_bucket_;

    /// Holds a bounding box for all possible shift vectors.
    PositionBoundingBoxType shift_bounding_box_;

    /// Holds an enlarged bounding box for all shift vectors.  It is
    /// larger by about half of a bucket in all directions.
    PositionBoundingBoxType shift_bounding_box_enlarged_;

    /// Diagonal size of each bucket in shift_bucket_.
    PositionType shift_bucket_size_;

    /// Number of surrounding buckets of element indices to be considered when
    /// computing shifts.
    Size element_bucket_window_[2];

    /// Number of surrounding buckets of shift indices to be considered when
    /// computing shifts.
    Size shift_bucket_window_[2];
  }
  ; // PoseClusteringShiftSuperimposer

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_PoseClusteringShiftSuperimposer_H
