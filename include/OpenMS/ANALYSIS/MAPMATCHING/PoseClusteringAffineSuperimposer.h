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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H

#include <OpenMS/KERNEL/DPeakConstReferenceArray.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>

#include <fstream>
#include <vector>
#include <map>
#include <math.h>

#define V_PoseClusteringAffineSuperimposer(bla) //  std::cout << bla << std::endl;

namespace OpenMS
{

  /**
    @brief Superimposer that uses a voting scheme to find a good affine transformation.

	  It works on two element maps (FeatureMap is the default map type, 
	  but you can also use a pointer map like DPeakConstReferenceArray) and 
	  computes a affine transformation, that maps the elements of one map (scene map) 
	  as near as possible to the elements in the other map (model map).
	  A element can be a DPeak, a DFeature, a ConsensusPeak or ConsensusFeature 
	  (wheras DFeature is the default element type).
	
	  This superimposer hashs all possible affine transformations and defines the 
	  transformation with the most votes as the best one.
		 
		@ref PoseClusteringAffineSuperimposer_Parameters are explained on a separate page.        
  */
  template < typename MapT = FeatureMap<> >
  class PoseClusteringAffineSuperimposer
        : public BaseSuperimposer< MapT >
  {
  public:

    enum HashMap
    {
      SHIFT = 0,
      SCALING = 1
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
  public:
    typedef typename Base::QualityType QualityType;
    typedef typename Base::PositionType PositionType;
    typedef typename Base::IntensityType IntensityType;
    typedef typename Base::PointType PointType;
    typedef typename Base::PointMapType PointMapType;
    typedef typename PositionType::CoordinateType CoordinateType;
    typedef DPeakConstReferenceArray< PointMapType > PeakPointerArray;
    typedef DBoundingBox<2>  PositionBoundingBoxType;
    typedef DBoundingBox<1> IntensityBoundingBoxType;
    typedef LinearMapping AffineTransformationType;
    typedef std::pair<int,int> PairType;
    typedef std::map< PairType, QualityType> AffineTransformationMapType;

    using Base::param_;
    using Base::setParameters;
    using Base::getParameters;
    using Base::subsections_;
    using Base::defaultsToParam_;
    using Base::defaults_;
    using Base::element_map_;
    using Base::final_transformation_;

    /// Constructor
    PoseClusteringAffineSuperimposer()
        : Base()
    {
      setName(getProductName());

      defaults_.setValue("tuple_search:mz_bucket_size",1,"An estimate of m/z deviation of corresponding elements in different maps.");
      defaults_.setValue("transformation_space:shift_bucket_size:RT",1,"Defines the shift parameter's bucket size during histograming.");
      defaults_.setValue("transformation_space:shift_bucket_size:MZ",0.1,"Defines the shift parameter's bucket size during histograming.");
      defaults_.setValue("transformation_space:scaling_bucket_size:RT",0.5,"Defines the scaling parameter's bucket size during histograming.");
      defaults_.setValue("transformation_space:scaling_bucket_size:MZ",0.1,"Defines the scaling parameter's bucket size during histograming.");
      defaults_.setValue("transformation_space:bucket_window_shift:RT",1,"Number of surrounding buckets of element indices to be considered when computing the shift parameter.");
      defaults_.setValue("transformation_space:bucket_window_shift:MZ",1,"Number of surrounding buckets of element indices to be considered when computing the shift parameter.");
      defaults_.setValue("transformation_space:bucket_window_scaling:RT",1,"Number of surrounding buckets of element indices to be considered when computing the scaling parameter.");
      defaults_.setValue("transformation_space:bucket_window_scaling:MZ",1,"Number of surrounding buckets of element indices to be considered when computing the scaling parameter.");
      defaults_.setValue("transformation_space:min_shift:RT",-1000,"Minimal shift parameter which is considered during histogramming.");
      defaults_.setValue("transformation_space:min_shift:MZ",-5,"Minimal shift parameter which is considered during histogramming.");
      defaults_.setValue("transformation_space:max_shift:RT",1000,"Maximal shift parameter which is considered during histogramming.");
      defaults_.setValue("transformation_space:max_shift:MZ",5,"Maximal shift parameter which is considered during histogramming.");
      defaults_.setValue("transformation_space:min_scaling:RT",-3,"Minimal scaling parameter which is considered during histogramming.");
      defaults_.setValue("transformation_space:min_scaling:MZ",-1.5,"Minimal scaling parameter which is considered during histogramming.");
      defaults_.setValue("transformation_space:max_scaling:RT",3,"Maximal scaling parameter which is considered during histogramming.");
      defaults_.setValue("transformation_space:max_scaling:MZ",1.5,"Maximal scaling parameter which is considered during histogramming.");

      defaultsToParam_();
    }

    /// Destructor
    virtual ~PoseClusteringAffineSuperimposer()
    {
      V_PoseClusteringAffineSuperimposer("~PoseClusteringAffineSuperimposer");
    }

    /// Estimates the transformation for each grid cell
    virtual void run()
    {
      if ( !this->element_map_[MODEL]->empty() && !this->element_map_[SCENE]->empty() )
      {
        // compute total intensities of both maps for normalisation
        IntensityType total_int_model_map = 0;
        UInt n = element_map_[MODEL]->size();
        for (UInt i = 0; i < n; ++i)
        {
          total_int_model_map += (*element_map_[MODEL])[i].getIntensity();
        }

        IntensityType total_int_scene_map = 0;
        n = element_map_[SCENE]->size();
        for (UInt i = 0; i < n; ++i)
        {
          total_int_scene_map += (*element_map_[SCENE])[i].getIntensity();
        }
        // clear scene_map_partners_
        scene_map_partners_.clear();
        // clear model_map_red_
        model_map_red_.clear();
        // clear rt hash
        rt_hash_.clear();
        // clear mz hash
        mz_hash_.clear();



        preprocess_();
        hashAffineTransformations_(total_int_model_map,total_int_scene_map);
        estimateFinalAffineTransformation_();
      }
      else
      {
        std::cerr << "PoseClusteringAffineSuperimposer::run():  Oops, one of the element maps is empty!\n";
      }
    }

    /// Returns an instance of this class
    static BaseSuperimposer<PointMapType>* create()
    {
      return new PoseClusteringAffineSuperimposer();
    }

    /// Returns the name of this module
    static const String getProductName()
    {
      return "poseclustering_affine";
    }

    /// Set size of the mz tolerance of point partners
    void setMzBucketSize(double mz_bucket_size)
    {
      mz_bucket_size_ = mz_bucket_size;
      param_.setValue("tuple_search:mz_bucket_size", mz_bucket_size);
    }

    /// Get size of the mz tolerance of point partners
    double getMzBucketSize() const
    {
      return mz_bucket_size_;
    }

    /// Set size of shift buckets (in dimension dim)
    void setShiftBucketSize(UInt dim, double shift_bucket_size)
    {
      shift_bucket_size_[dim] = shift_bucket_size;
      param_.setValue(String("transformation_space:shift_bucket_size:") + RawDataPoint2D::shortDimensionName(dim), shift_bucket_size);
    }

    /// Get size of shift buckets (in dimension dim)
    double getShiftBucketSize(UInt dim) const
    {
      return shift_bucket_size_[dim];
    }

    /// Set size of scaling buckets (in dimension dim)
    void setScalingBucketSize(UInt dim, double scaling_bucket_size)
    {
      scaling_bucket_size_[dim] = scaling_bucket_size;
      param_.setValue(String("transformation_space:scaling_bucket_size:") + RawDataPoint2D::shortDimensionName(dim), scaling_bucket_size);
    }

    /// Get size of scaling buckets (in dimension dim)
    double getScalingBucketSize(UInt dim) const
    {
      return scaling_bucket_size_[dim];
    }

    /// Set number of neighbouring shift buckets to be considered for the calculation of the final transformation (in dimension dim)
    void setBucketWindowShift(UInt dim, UInt bucket_window_shift)
    {
      bucket_window_shift_[dim] = bucket_window_shift;
      param_.setValue(String("transformation_space:bucket_window_shift:") + RawDataPoint2D::shortDimensionName(dim), (int)bucket_window_shift);
    }

    /// Get number of neighbouring shift buckets to be considered for the calculation of the final transformation (in dimension dim)
    UInt getBucketWindowShift(UInt dim) const
    {
      return bucket_window_shift_[dim];
    }

    /// Set number of neighbouring scaling buckets to be considered for the calculation of the final transformation (in dimension dim)
    void setBucketWindowScaling(UInt dim, UInt bucket_window_scaling)
    {
      bucket_window_scaling_[dim] = bucket_window_scaling;
      param_.setValue(String("transformation_space:bucket_window_scaling:") + RawDataPoint2D::shortDimensionName(dim), (int)bucket_window_scaling);
    }

    /// Get number of neighbouring scaling buckets to be considered for the calculation of the final transformation (in dimension dim)
    UInt getBucketWindowScaling(UInt dim) const
    {
      return bucket_window_scaling_[dim];
    }


  protected:
    virtual void updateMembers_()
    {
      mz_bucket_size_ = (CoordinateType)param_.getValue("tuple_search:mz_bucket_size");
      shift_bucket_size_[0] = (CoordinateType)param_.getValue("transformation_space:shift_bucket_size:RT");
      shift_bucket_size_[1] = (CoordinateType)param_.getValue("transformation_space:shift_bucket_size:MZ");
      scaling_bucket_size_[0] = (CoordinateType)param_.getValue("transformation_space:scaling_bucket_size:RT");
      scaling_bucket_size_[1] = (CoordinateType)param_.getValue("transformation_space:scaling_bucket_size:MZ");
      bucket_window_shift_[0]  = (UInt)param_.getValue("transformation_space:bucket_window_shift:RT");
      bucket_window_shift_[1]  = (UInt)param_.getValue("transformation_space:bucket_window_shift:MZ");
      bucket_window_scaling_[0] = (UInt)param_.getValue("transformation_space:bucket_window_scaling:RT");
      bucket_window_scaling_[1] = (UInt)param_.getValue("transformation_space:bucket_window_scaling:MZ");

      PositionType min;
      PositionType max;
      min[RawDataPoint2D::RT] = (CoordinateType)param_.getValue("transformation_space:min_shift:RT");
      min[RawDataPoint2D::MZ] = (CoordinateType)param_.getValue("transformation_space:min_shift:MZ");
      max[RawDataPoint2D::RT] = (CoordinateType)param_.getValue("transformation_space:max_shift:RT");
      max[RawDataPoint2D::MZ] = (CoordinateType)param_.getValue("transformation_space:max_shift:MZ");

      shift_bounding_box_.enlarge(min);
      shift_bounding_box_.enlarge(max);

      min[RawDataPoint2D::RT] = (CoordinateType)param_.getValue("transformation_space:min_scaling:RT");
      min[RawDataPoint2D::MZ] = (CoordinateType)param_.getValue("transformation_space:min_scaling:MZ");
      max[RawDataPoint2D::RT] = (CoordinateType)param_.getValue("transformation_space:max_scaling:RT");
      max[RawDataPoint2D::MZ] = (CoordinateType)param_.getValue("transformation_space:max_scaling:MZ");

      scaling_bounding_box_.enlarge(min);
      scaling_bounding_box_.enlarge(max);

    }

    /// To speed up the calculation of the final transformation, we confine the number of
    /// considered point pairs. We match a point p in the model map only onto those points p'
    /// in the scene map that lie in a certain mz intervall.
    /// If  (p_mz - mz_bucket_size_) <= p'_mz <= (p_mz mz_bucket_size_) then p and p' are partners.
    void preprocess_()
    {
#define V_preprocessSceneMap(bla)  V_PoseClusteringAffineSuperimposer(bla)
      // build an array of pointer to all elements in the model map
      PeakPointerArray model_map(element_map_[MODEL]->begin(), element_map_[MODEL]->end());
      // build an array of pointer to all elements in the scene map
      PeakPointerArray scene_map(element_map_[SCENE]->begin(), element_map_[SCENE]->end());

      // for each element (rt_m,mz_m) of the model map
      // search for corresponding elements (rt_i,mz_i) in the scene map
      // which lie in a predefined mz intervall (mz_i in [mz_m-eps,mz_m+eps))
      model_map.sortByNthPosition(RawDataPoint2D::MZ);
      scene_map.sortByNthPosition(RawDataPoint2D::MZ);

      // take only elements of the model map which have partners in the scene map
      typename PeakPointerArray::const_iterator it_first = scene_map.begin();
      typename PeakPointerArray::const_iterator it_last = it_first;
      UInt n = model_map.size();
      for (UInt i = 0; i < n; ++i)
      {
        DoubleReal act_mz = model_map[i].getMZ();
        DoubleReal min_mz = act_mz - mz_bucket_size_;
        DoubleReal max_mz = act_mz + mz_bucket_size_;

        CoordinateType act_rt = model_map[i].getRT();
        CoordinateType min_rt = act_rt - 5000;
        CoordinateType max_rt = act_rt + 5000;

        std::vector< const PointType* > partners;
        // search for the left end of the intervall
        while ((it_first >= scene_map.begin()) && (it_first != scene_map.end()) && (it_first->getMZ() < min_mz)
               && (it_first->getMZ() < max_mz) && (it_first->getMZ() < min_mz) )
        {
          ++it_first;
        }

        it_last = it_first;

        // search for the right end of the intervall
        while ((it_last < scene_map.end()) && (it_last->getMZ() < max_mz))
        {
          ++it_last;
        }

        /// TODO Remove
        PeakPointerArray partner_;
        for (typename PeakPointerArray::const_iterator it = it_first; it != it_last; ++it)
        {
          if ((it->getRT() < max_rt) && (it->getRT() > min_rt) )
          {
            partners.push_back(&(*it));
          }

        }

        if (partners.size() > 0)
        {
          model_map_red_.push_back(model_map[i]);
          scene_map_partners_.push_back(partners);
        }
      }
      //       std::cout.flush();

      // Compute shift_bucket_size_ and num_buckets.
      PositionType diagonal_shift = shift_bounding_box_.diagonal();
      PositionType diagonal_scaling = scaling_bounding_box_.diagonal();
      V_preprocessSceneMap("diagonal shift: " << diagonal_shift);
      V_preprocessSceneMap("diagonal scaling: " << diagonal_scaling);

      for ( UInt dimension = 0; dimension < 2; ++dimension)
      {
        num_buckets_shift_[dimension] = (int)(ceil(diagonal_shift[dimension]/shift_bucket_size_[dimension]));
        num_buckets_scaling_[dimension] = (int)(ceil(diagonal_scaling[dimension]/scaling_bucket_size_[dimension]));
        shift_bucket_size_[dimension] = diagonal_shift[dimension] / (CoordinateType)num_buckets_shift_[dimension];
        scaling_bucket_size_[dimension] = diagonal_scaling[dimension] / (CoordinateType)num_buckets_scaling_[dimension];
      }
      V_preprocessSceneMap("shift bucket size : " << shift_bucket_size_[RawDataPoint2D::RT] << ' ' << shift_bucket_size_[RawDataPoint2D::MZ]);
      V_preprocessSceneMap("scaling bucket size : " << scaling_bucket_size_);
      V_preprocessSceneMap("shift number_buckets: " << num_buckets_shift_[0] << ' ' << num_buckets_shift_[1]);
      V_preprocessSceneMap("scaling number_buckets : " << num_buckets_scaling_[0] << ' ' << num_buckets_scaling_[1]);
#undef V_preprocessSceneMap

    } // preprocess_

    /// Compute the transformations between each point pair in the model map and each point pair in the scene map
    /// and hash the affine transformation.
    void hashAffineTransformations_(IntensityType total_int_model_map, IntensityType total_int_scene_map )
    {
#define V_hashAffineTransformations_(bla) V_PoseClusteringAffineSuperimposer(bla)
      // take each point pair in the model map
      UInt n = model_map_red_.size();
      for (UInt i = 0; i < n; ++i)
      {
        // take only the next 10 neighbours in m/z as partner in the model map
        UInt k=((i+50)>= n) ? n : (i+50);
        for (UInt j = i+1; j < k; ++j)
        {
          // avoid cross mappings (i,j) -> (k,l) (e.g. i_rt < j_rt and k_rt > l_rt)
          // and point pairs with equal retention times (e.g. i_rt == j_rt)
          PositionType diff = model_map_red_[i].getPosition() - model_map_red_[j].getPosition();
          // and compute the affine transformation to all corresponding points pair in the scene map
          std::vector< const PointType* >& partners_i = scene_map_partners_[i];
          std::vector< const PointType* >& partners_j  = scene_map_partners_[j];
          UInt m = partners_i.size();
          UInt p = partners_j.size();

          for (UInt k = 0; k < m; ++k)
          {
            for (UInt l = 0; l < p; ++l)
            {
              //               std::cout << "Partner of " << model_map_red_[j].getPosition() << ' ' << l << ' ' << partners_j.size() << std::endl;
              //               std::cout << *(partners_j[l]) << std::endl;
              PositionType diff_2 = (partners_j[l]->getPosition() - partners_i[k]->getPosition());

              // compute the transformation (i,j) -> (k,l)
              PositionType shift;
              PositionType scaling;
              bool transformation_ok[2];
              transformation_ok[0] = false;
              transformation_ok[1] = false;

              if ((fabs(diff[RawDataPoint2D::RT]) > 0.001) && (fabs(diff_2[RawDataPoint2D::RT]) > 0.001))
              {
                scaling[RawDataPoint2D::RT] = (model_map_red_[j].getRT() - model_map_red_[i].getRT())
                                              /(partners_j[l]->getRT() - partners_i[k]->getRT());

                shift[RawDataPoint2D::RT] =  model_map_red_[i].getRT() - partners_i[k]->getRT()*scaling[RawDataPoint2D::RT];
                transformation_ok[RawDataPoint2D::RT] = true;
              }


              if ((fabs(diff[RawDataPoint2D::MZ]) > 0.001) && (fabs(diff_2[RawDataPoint2D::MZ]) > 0.001))
              {
                scaling[RawDataPoint2D::MZ] = (model_map_red_[j].getMZ() - model_map_red_[i].getMZ())
                                              /(partners_j[l]->getMZ() - partners_i[k]->getMZ());

                shift[RawDataPoint2D::MZ] =  model_map_red_[i].getMZ() - partners_i[k]->getMZ()*scaling[RawDataPoint2D::MZ];
                transformation_ok[RawDataPoint2D::MZ] = true;
              }

              // compute the hash indices
              int bucket_shift_index;
              int bucket_scaling_index;
              CoordinateType bucket_fraction_shift = 0.;
              CoordinateType bucket_fraction_scaling = 0.;

              IntensityType int_i = model_map_red_[i].getIntensity()/total_int_model_map;
              IntensityType int_k = partners_i[k]->getIntensity()/total_int_scene_map;
              IntensityType int_j = model_map_red_[j].getIntensity()/total_int_model_map;
              IntensityType int_l = partners_j[l]->getIntensity()/total_int_scene_map;
              QualityType int_similarity_jl = (int_j < int_l) ? int_j/int_l : int_l/int_j;
              QualityType int_similarity_ik = (int_i < int_k) ? int_i/int_k : int_k/int_i;
              QualityType int_similarity =  int_similarity_ik + int_similarity_jl;

              // check if the transformation paramerters lie in the hash map
              if ( transformation_ok[RawDataPoint2D::RT] && (scaling_bounding_box_.min()[RawDataPoint2D::RT] < scaling[RawDataPoint2D::RT]) && (scaling_bounding_box_.max()[RawDataPoint2D::RT] > scaling[RawDataPoint2D::RT])
                   && (shift_bounding_box_.min()[RawDataPoint2D::RT] < shift[RawDataPoint2D::RT]) && (shift_bounding_box_.max()[RawDataPoint2D::RT] > shift[RawDataPoint2D::RT]))
              {
                bucket_fraction_shift = (shift[RawDataPoint2D::RT] - shift_bounding_box_.min()[RawDataPoint2D::RT]) / shift_bucket_size_[RawDataPoint2D::RT];  // floating point division
                bucket_fraction_scaling = (scaling[RawDataPoint2D::RT] - scaling_bounding_box_.min()[RawDataPoint2D::RT]) / scaling_bucket_size_[RawDataPoint2D::RT];  // floating point division
                bucket_shift_index    = (int) bucket_fraction_shift; // round down (yes we are >= 0)
                //                       std::cout << "bucket_shift_index " <<  bucket_shift_index[RawDataPoint2D::RT] << std::endl;
                bucket_scaling_index    = (int) bucket_fraction_scaling; // round down (yes we are >= 0)
                //                       std::cout << "bucket_scaling_index " << bucket_shift_index[RawDataPoint2D::RT] << std::endl;
                bucket_fraction_shift -= bucket_shift_index;          // fractional part
                //                       std::cout << "bucket_fraction_shift " << bucket_fraction_shift[RawDataPoint2D::RT]<< std::endl;
                bucket_fraction_scaling -= bucket_scaling_index;          // fractional part
                //                       std::cout << "bucket_fraction_scaling" << bucket_fraction_scaling[RawDataPoint2D::RT]<< std::endl;

                //                     std::cout << "Buckets rt (" << bucket_shift_index[RawDataPoint2D::RT] << ',' << bucket_scaling_index[RawDataPoint2D::RT] << ')' << std::endl;
                //                     std::cout << "Buckets mz (" << bucket_shift_index[RawDataPoint2D::MZ] << ',' << bucket_scaling_index[RawDataPoint2D::MZ] << ')' << std::endl;

                // hash the transformation if possible
                CoordinateType bucket_fraction_complement_shift = 1.;
                CoordinateType bucket_fraction_complement_scaling = 1.;
                bucket_fraction_complement_shift -= bucket_fraction_shift;
                bucket_fraction_complement_scaling -= bucket_fraction_scaling;

                // Distribute the vote of the shift among the four neighboring buckets.
                QualityType factor;
                // rt
                factor = bucket_fraction_complement_shift * bucket_fraction_complement_scaling * int_similarity;
                rt_hash_[PairType(bucket_shift_index, bucket_scaling_index)] += factor;
                //                     std::cout << "hash shift " << bucket_shift_index*shift_bucket_size_ + shift_bounding_box_.min() << std::endl;
                //                     std::cout << "hash scaling " << bucket_scaling_index*scaling_bucket_size_ + scaling_bounding_box_.min()[RawDataPoint2D::RT] << std::endl;
                //                     std::cout << "factor" << factor << std::endl;

                factor = bucket_fraction_complement_shift * bucket_fraction_scaling * int_similarity;
                rt_hash_[PairType(bucket_shift_index, bucket_scaling_index + 1 )] += factor;
                //                     std::cout << "hash shift " << bucket_shift_index*shift_bucket_size_ + shift_bounding_box_.min()[RawDataPoint2D::RT] << std::endl;
                //                     std::cout << "hash scaling " << (bucket_scaling_index+1)*scaling_bucket_size_ + scaling_bounding_box_.min() << std::endl;

                factor = bucket_fraction_shift * bucket_fraction_complement_scaling * int_similarity;
                rt_hash_[PairType( bucket_shift_index + 1, bucket_scaling_index )] += factor;

                factor = bucket_fraction_shift * bucket_fraction_scaling * int_similarity_ik;
                rt_hash_[PairType( bucket_shift_index + 1, bucket_scaling_index + 1 )] += factor;
              }

              if ( transformation_ok[RawDataPoint2D::MZ] && (scaling_bounding_box_.min()[RawDataPoint2D::MZ] < scaling[RawDataPoint2D::MZ]) && (scaling_bounding_box_.max()[RawDataPoint2D::MZ] > scaling[RawDataPoint2D::MZ])
                   && (shift_bounding_box_.min()[RawDataPoint2D::MZ] < shift[RawDataPoint2D::MZ]) && (shift_bounding_box_.max()[RawDataPoint2D::MZ] > shift[RawDataPoint2D::MZ]))
              {
                bucket_fraction_shift = (shift[RawDataPoint2D::MZ] - shift_bounding_box_.min()[RawDataPoint2D::MZ]) / shift_bucket_size_[RawDataPoint2D::MZ];  // floating point division
                bucket_fraction_scaling = (scaling[RawDataPoint2D::MZ] - scaling_bounding_box_.min()[RawDataPoint2D::MZ]) / scaling_bucket_size_[RawDataPoint2D::MZ];  // floating point division
                bucket_shift_index    = (int) bucket_fraction_shift; // round down (yes we are >= 0)
                //                       std::cout << "bucket_shift_index " <<  bucket_shift_index[RawDataPoint2D::MZ] << std::endl;
                bucket_scaling_index    = (int) bucket_fraction_scaling; // round down (yes we are >= 0)
                //                       std::cout << "bucket_scaling_index " << bucket_shift_index[RawDataPoint2D::MZ] << std::endl;
                bucket_fraction_shift -= bucket_shift_index;          // fractional part
                //                       std::cout << "bucket_fraction_shift " << bucket_fraction_shift[RawDataPoint2D::MZ]<< std::endl;
                bucket_fraction_scaling -= bucket_scaling_index;          // fractional part
                //                       std::cout << "bucket_fraction_scaling" << bucket_fraction_scaling[RawDataPoint2D::MZ]<< std::endl;

                //                     std::cout << "Buckets rt (" << bucket_shift_index[RawDataPoint2D::MZ] << ',' << bucket_scaling_index[RawDataPoint2D::MZ] << ')' << std::endl;
                //                     std::cout << "Buckets mz (" << bucket_shift_index[RawDataPoint2D::MZ] << ',' << bucket_scaling_index[RawDataPoint2D::MZ] << ')' << std::endl;

                // hash the transformation if possible
                CoordinateType bucket_fraction_complement_shift = 1.;
                CoordinateType bucket_fraction_complement_scaling = 1.;
                bucket_fraction_complement_shift -= bucket_fraction_shift;
                bucket_fraction_complement_scaling -= bucket_fraction_scaling;

                // Distribute the vote of the shift among the four neighboring buckets.
                QualityType factor;

                // mz
                factor = bucket_fraction_complement_shift * bucket_fraction_complement_scaling * int_similarity;
                mz_hash_[PairType(bucket_shift_index, bucket_scaling_index)] += factor;
                //                     std::cout << "hash shift " << bucket_shift_index*shift_bucket_size_ + shift_bounding_box_.min() << std::endl;
                //                     std::cout << "hash scaling " << bucket_scaling_index*scaling_bucket_size_ + scaling_bounding_box_.min()[RawDataPoint2D::MZ] << std::endl;
                //                     std::cout << "factor" << factor << std::endl;

                factor = bucket_fraction_complement_shift * bucket_fraction_scaling * int_similarity;
                mz_hash_[PairType(bucket_shift_index, bucket_scaling_index + 1 )] += factor;
                //                     std::cout << "hash shift " << bucket_shift_index*shift_bucket_size_ + shift_bounding_box_.min()[RawDataPoint2D::MZ] << std::endl;
                //                     std::cout << "hash scaling " << (bucket_scaling_index+1)*scaling_bucket_size_ + scaling_bounding_box_.min() << std::endl;

                factor = bucket_fraction_shift * bucket_fraction_complement_scaling * int_similarity;
                mz_hash_[PairType( bucket_shift_index + 1, bucket_scaling_index )] += factor;

                factor = bucket_fraction_shift * bucket_fraction_scaling * int_similarity;
                mz_hash_[PairType( bucket_shift_index + 1, bucket_scaling_index + 1 )] += factor;
              }
            } // for l
          } // for k
        } // for j
      } // for i



      //       std::ofstream rt_os("rt_matrix.dat", std::ios::out);
      //       std::ofstream mz_os("mz_matrix.dat", std::ios::out);
      //
      //       typename AffineTransformationMapType::const_iterator it = rt_hash_.begin();
      //       while (it != rt_hash_.end())
      //       {
      //         rt_os << ((it->first).first)*shift_bucket_size_[RawDataPoint2D::RT] + shift_bounding_box_.min()[RawDataPoint2D::RT] << ' '
      //         << ((it->first).second)*scaling_bucket_size_[RawDataPoint2D::RT] + scaling_bounding_box_.min()[RawDataPoint2D::RT] << ' '
      //         << it->second << '\n';
      //         ++it;
      //       }
      //
      //       it = mz_hash_.begin();
      //       while (it != mz_hash_.end())
      //       {
      //         mz_os << ((it->first).first)*shift_bucket_size_[RawDataPoint2D::MZ] + shift_bounding_box_.min()[RawDataPoint2D::MZ] << ' '
      //         << ((it->first).second)*scaling_bucket_size_[RawDataPoint2D::MZ] + scaling_bounding_box_.min()[RawDataPoint2D::MZ] << ' '
      //         << it->second << '\n';
      //         ++it;
      //       }
      //       rt_os.flush();
      //       mz_os.flush();
#undef V_hashAffineTransformations_

    } // hashAffineTransformations_


    /// After the hashing phase, the best transformation, that is the transformation with the most votes
    /// is determined.
    void estimateFinalAffineTransformation_()
    {
#define V_estimateFinalAffineTransformation_(bla) V_PoseClusteringAffineSuperimposer(bla)
      V_estimateFinalAffineTransformation_("@@@ computeShift_()");

      // search for the maximal vote parameter of the rt transformation
      PairType max_element_index_rt;
      QualityType act_max_rt = 0;
      for (typename AffineTransformationMapType::const_iterator it = rt_hash_.begin(); it != rt_hash_.end(); ++it)
      {
        if (it->second > act_max_rt)
        {
          max_element_index_rt = it->first;
          act_max_rt = it->second;
        }
      }

      V_estimateFinalAffineTransformation_("Max element in rt: Indizes: "<< max_element_index_rt.first << ' ' << max_element_index_rt.second
                                           << " Votes: " << act_max_rt
                                           << " shift: "  << max_element_index_rt.first*shift_bucket_size_[RawDataPoint2D::RT] + shift_bounding_box_.min()[RawDataPoint2D::RT]
                                           << " scaling: " << max_element_index_rt.second*scaling_bucket_size_[RawDataPoint2D::RT] + scaling_bounding_box_.min()[RawDataPoint2D::RT]);


      // Compute a weighted average of the transformation parameters nearby the max_element_index.
      PositionType rt_trafo;
      PositionType rt_bounding_box_min(shift_bounding_box_.min()[RawDataPoint2D::RT],scaling_bounding_box_.min()[RawDataPoint2D::RT]);
      int rt_run_indices[2];
      QualityType quality = 0;
      PositionType rt_window(bucket_window_shift_[RawDataPoint2D::RT],bucket_window_scaling_[RawDataPoint2D::RT]);
      //         std::cout << "rt_window " << rt_window << std::endl;
      for ( rt_run_indices[SHIFT]  = std::max ( int (max_element_index_rt.first - rt_window[SHIFT]), 0 );
            rt_run_indices[SHIFT] <= std::min ( int (max_element_index_rt.first + rt_window[SHIFT]), num_buckets_shift_[RawDataPoint2D::RT] );
            ++rt_run_indices[SHIFT])
      {
        for ( rt_run_indices[SCALING]  = std::max ( int (max_element_index_rt.second - rt_window[SCALING]), 0 );
              rt_run_indices[SCALING] <= std::min ( int (max_element_index_rt.second + rt_window[SCALING]), num_buckets_scaling_[RawDataPoint2D::RT] );
              ++rt_run_indices[SCALING])

        {
          PositionType contribution_position(shift_bucket_size_[RawDataPoint2D::RT],scaling_bucket_size_[RawDataPoint2D::RT]);
          // is the neighbouring bucket in the map?
          typename AffineTransformationMapType::const_iterator it = rt_hash_.find(PairType(rt_run_indices[SHIFT], rt_run_indices[SCALING]));
          if ( it != rt_hash_.end())
          {
            for ( UInt dimension = 0; dimension < 2; ++dimension)
            {
              contribution_position[dimension] *= rt_run_indices[dimension];
            }
            contribution_position += rt_bounding_box_min;
            // std::cout << "contribution_position " << contribution_position << std::endl;
            QualityType contribution_quality = it->second;
            // std::cout << "contribution_quality " << contribution_quality << std::endl;
            quality += contribution_quality;
            // std::cout << "quality " << quality << std::endl;
            contribution_position *= contribution_quality;
            // std::cout << "contribution_position " << contribution_position << std::endl;
            rt_trafo += contribution_position;
            // std::cout << "rt_trafo " << rt_trafo << std::endl;
          }
        }
      }

      if ( quality != 0 )
      {
        rt_trafo /= quality;
      }

      // Assign the result.
      // set slope and intercept
      final_transformation_[RawDataPoint2D::RT].setParam( rt_trafo[SCALING], rt_trafo[SHIFT] );
      V_estimateFinalAffineTransformation_("estimateFinalAffineTransformation_() hat geklappt rt: " << rt_trafo);


      // search for the maximal vote parameter of the mz transformation
      PairType max_element_index_mz;
      QualityType act_max_mz = 0;
      for (typename AffineTransformationMapType::const_iterator it = mz_hash_.begin(); it != mz_hash_.end(); ++it)
      {
        if (it->second > act_max_mz)
        {
          max_element_index_mz = it->first;
          act_max_mz = it->second;
        }
      }

      V_estimateFinalAffineTransformation_("Max element in mz: Indizes: "<< max_element_index_mz.first << ' ' << max_element_index_mz.second
                                           << " Votes: " << act_max_mz
                                           << " shift: "  << max_element_index_mz.first*shift_bucket_size_[RawDataPoint2D::MZ] + shift_bounding_box_.min()[RawDataPoint2D::MZ]
                                           << " scaling: " << max_element_index_mz.second*scaling_bucket_size_[RawDataPoint2D::MZ] + scaling_bounding_box_.min()[RawDataPoint2D::MZ]);

      PositionType mz_trafo;
      PositionType mz_bounding_box_min(shift_bounding_box_.min()[RawDataPoint2D::MZ],scaling_bounding_box_.min()[RawDataPoint2D::MZ]);
      int mz_run_indices[2];
      quality=0;
      PositionType mz_window(bucket_window_shift_[RawDataPoint2D::MZ],bucket_window_scaling_[RawDataPoint2D::MZ]);
      for ( mz_run_indices[SHIFT]  = std::max ( int (max_element_index_mz.first - mz_window[SHIFT]), 0 );
            mz_run_indices[SHIFT] <= std::min ( int (max_element_index_mz.first + mz_window[SHIFT]), num_buckets_shift_[RawDataPoint2D::MZ] );
            ++mz_run_indices[SHIFT])
      {
        for ( mz_run_indices[SCALING]  = std::max ( int (max_element_index_mz.second - mz_window[SCALING]), 0 );
              mz_run_indices[SCALING] <= std::min ( int (max_element_index_mz.second + mz_window[SCALING]), num_buckets_scaling_[RawDataPoint2D::MZ] );
              ++mz_run_indices[SCALING])
        {
          PositionType contribution_position(shift_bucket_size_[RawDataPoint2D::MZ],scaling_bucket_size_[RawDataPoint2D::MZ]);
          // is the neighbouring bucket in the map?
          typename AffineTransformationMapType::const_iterator it = mz_hash_.find(PairType(mz_run_indices[SHIFT], mz_run_indices[SCALING]));
          if ( it != mz_hash_.end())
          {
            for ( UInt dimension = 0; dimension < 2; ++dimension)
            {
              contribution_position[dimension] *= mz_run_indices[dimension];
            }
            contribution_position += mz_bounding_box_min;
            //               std::cout << "contribution_position " << contribution_position << std::endl;
            QualityType contribution_quality = it->second;
            //               std::cout << "contribution_quality " << contribution_quality << std::endl;
            quality += contribution_quality;
            //               std::cout << "quality " << quality << std::endl;
            contribution_position *= contribution_quality;
            //               std::cout << "contribution_position " << contribution_position << std::endl;
            mz_trafo += contribution_position;
            //               std::cout << "mz_trafo " << mz_trafo << std::endl;
          }
        }
      }
      if ( quality != 0 )
      {
        mz_trafo /= quality;
      }
      // set slope and intercept
      final_transformation_[RawDataPoint2D::MZ].setParam( mz_trafo[SCALING], mz_trafo[SHIFT] );
      V_estimateFinalAffineTransformation_("estimateFinalAffineTransformation_() hat geklappt mz: " << mz_trafo);

      //         std::ofstream rt_os("rt_matrix.pgm", std::ios::out);
      //         std::ofstream mz_os("mz_matrix.pgm", std::ios::out);
      //
      //         typename AffineTransformationMapType::const_iterator it = rt_hash_.begin();
      //         while (it != rt_hash_.end())
      //         {
      //           rt_os << (it->first).first*shift_bucket_size_[RawDataPoint2D::RT] + shift_bounding_box_.min()[RawDataPoint2D::RT]
      //           << ' ' << (it->first).second*scaling_bucket_size_[RawDataPoint2D::RT] + scaling_bounding_box_.min()[RawDataPoint2D::RT] << '\n';
      //           ++it;
      //         }
      //         rt_os.flush();
#undef V_estimateFinalAffineTransformation_

    } // estimateFinalAffineTransformation_

    /// Reduced model map which contains only elements of the model map which have a partner in the scene map
    PeakPointerArray model_map_red_;

    /// Partner elements in the scene map
    std::vector< std::vector< const PointType* > > scene_map_partners_;

    /// Hash map of all transformations in the rt dimension (contains the affine transformation parameters)
    AffineTransformationMapType rt_hash_;

    /// Hash map of all transformations in the mz dimension (contains the affine transformation parameters)
    AffineTransformationMapType mz_hash_;

    /// Bounding box of the shift parameters
    PositionBoundingBoxType shift_bounding_box_;

    /// Bounding box of the scaling parameters
    PositionBoundingBoxType scaling_bounding_box_;

    /// Diagonal size of each shift bucket
    PositionType shift_bucket_size_;

    /// Diagonal size of each scaling bucket
    PositionType scaling_bucket_size_;

    /// Number of shift buckets
    int num_buckets_shift_[2];

    /// Number of scaling buckets
    int num_buckets_scaling_[2];

    /// Number of neighbouring shift buckets to be considered computing the final transformations
    UInt bucket_window_shift_[2];

    /// Number of neighbouring scaling buckets to be considered computing the final transformations
    UInt bucket_window_scaling_[2];

    /// Maximum deviation in mz of two partner points
    CoordinateType mz_bucket_size_;

  }
  ; // PoseClusteringAffineSuperimposer
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_PoseClusteringAffineSuperimposer_H
