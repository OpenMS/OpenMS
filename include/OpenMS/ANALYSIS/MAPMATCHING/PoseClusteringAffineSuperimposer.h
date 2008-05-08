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


#ifndef OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H

#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/DPeakConstReferenceArray.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/LinearMapping.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>

#include <fstream>
#include <vector>
#include <map>
#include <math.h>

namespace OpenMS
{

  /**
    @brief Superimposer that uses a voting scheme to find a good affine transformation.

	  It works on two element maps (FeatureMap is the default map type) and 
	  computes a affine transformation, that maps the elements of one map (scene map) 
	  as near as possible to the elements in the other map (model map).
	  A element can be a DPeak, a DFeature or ConsensusFeature 
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
  		
  		///Base class type definition
			typedef BaseSuperimposer< MapT > Base;

	    typedef typename Base::ElementMapType ElementMapType;
	    typedef typename ElementMapType::value_type PointType;
	    typedef DPeakConstReferenceArray<ElementMapType> PeakPointerArray;
	    typedef std::pair<int,int> PairType;
	    typedef std::map< PairType, DoubleReal> AffineTransformationMapType;
	
	    using Base::param_;
	    using Base::defaultsToParam_;
	    using Base::defaults_;
	    using Base::model_map_;
	    using Base::scene_map_;

	    /// Constructor
	    PoseClusteringAffineSuperimposer()
	    	: Base()
	    {
	      Base::setName(getProductName());
	
	      defaults_.setValue("tuple_search:mz_bucket_size",.5,"An estimate of m/z deviation of corresponding elements in different maps.");
	      defaults_.setValue("transformation_space:shift_bucket_size:RT",10.0,"Defines the shift parameter's bucket size during histograming.",true);
	      defaults_.setValue("transformation_space:shift_bucket_size:MZ",0.01,"Defines the shift parameter's bucket size during histograming.",true);
	      defaults_.setValue("transformation_space:scaling_bucket_size:RT",0.01,"Defines the scaling parameter's bucket size during histograming.",true);
	      defaults_.setValue("transformation_space:scaling_bucket_size:MZ",0.01,"Defines the scaling parameter's bucket size during histograming.",true);
	      defaults_.setValue("transformation_space:bucket_window_shift:RT",2,"Number of surrounding buckets of element indices to be considered when computing the shift parameter.",true);
	      defaults_.setValue("transformation_space:bucket_window_shift:MZ",2,"Number of surrounding buckets of element indices to be considered when computing the shift parameter.",true);
	      defaults_.setValue("transformation_space:bucket_window_scaling:RT",2,"Number of surrounding buckets of element indices to be considered when computing the scaling parameter.",true);
	      defaults_.setValue("transformation_space:bucket_window_scaling:MZ",2,"Number of surrounding buckets of element indices to be considered when computing the scaling parameter.",true);
	      defaults_.setValue("transformation_space:min_shift:RT",-1000.0,"Minimal shift parameter which is considered during histogramming.",true);
	      defaults_.setValue("transformation_space:min_shift:MZ",-5.0,"Minimal shift parameter which is considered during histogramming.",true);
	      defaults_.setValue("transformation_space:max_shift:RT",1000.0,"Maximal shift parameter which is considered during histogramming.",true);
	      defaults_.setValue("transformation_space:max_shift:MZ",5.0,"Maximal shift parameter which is considered during histogramming.",true);
	      defaults_.setValue("transformation_space:min_scaling:RT",0.5,"Minimal scaling parameter which is considered during histogramming.",true);
	      defaults_.setValue("transformation_space:min_scaling:MZ",0.8,"Minimal scaling parameter which is considered during histogramming.",true);
	      defaults_.setValue("transformation_space:max_scaling:RT",2.0,"Maximal scaling parameter which is considered during histogramming.",true);
	      defaults_.setValue("transformation_space:max_scaling:MZ",1.2,"Maximal scaling parameter which is considered during histogramming.",true);
	
	      defaultsToParam_();
	    }

	    /// Destructor
	    virtual ~PoseClusteringAffineSuperimposer()
	    {
	    }

	    /// Estimates the transformation for each grid cell
	    virtual void run(LinearMapping& mapping)
	    {
				if (model_map_==0) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"model_map");
				if (scene_map_==0) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"scene_map");
				
        // compute total intensities of both maps for normalisation
        DoubleReal total_int_model_map = 0;
        UInt n = model_map_->size();
        for (UInt i = 0; i < n; ++i)
        {
          total_int_model_map += (*model_map_)[i].getIntensity();
        }

        DoubleReal total_int_scene_map = 0;
        n = scene_map_->size();
        for (UInt i = 0; i < n; ++i)
        {
          total_int_scene_map += (*scene_map_)[i].getIntensity();
        }
        
        // clear 
        scene_map_partners_.clear();
        model_map_red_.clear();
        rt_hash_.clear();
        mz_hash_.clear();

        preprocess_();
        hashAffineTransformations_(total_int_model_map,total_int_scene_map);
        estimateFinalAffineTransformation_(mapping);
    	}

	    /// Returns an instance of this class
	    static BaseSuperimposer<ElementMapType>* create()
	    {
	      return new PoseClusteringAffineSuperimposer();
	    }
	
	    /// Returns the name of this module
	    static const String getProductName()
	    {
	      return "poseclustering_affine";
	    }
	
	  protected:
	  	
	    virtual void updateMembers_()
	    {
	      mz_bucket_size_ = (DoubleReal)param_.getValue("tuple_search:mz_bucket_size");
	      shift_bucket_size_ = (DoubleReal)param_.getValue("transformation_space:shift_bucket_size:RT");
	      scaling_bucket_size_ = (DoubleReal)param_.getValue("transformation_space:scaling_bucket_size:RT");
	      bucket_window_shift_  = (UInt)param_.getValue("transformation_space:bucket_window_shift:RT");
	      bucket_window_scaling_ = (UInt)param_.getValue("transformation_space:bucket_window_scaling:RT");
	
	      DPosition<2> min;
	      DPosition<2> max;
	      min[RawDataPoint2D::RT] = (DoubleReal)param_.getValue("transformation_space:min_shift:RT");
	      min[RawDataPoint2D::MZ] = (DoubleReal)param_.getValue("transformation_space:min_shift:MZ");
	      max[RawDataPoint2D::RT] = (DoubleReal)param_.getValue("transformation_space:max_shift:RT");
	      max[RawDataPoint2D::MZ] = (DoubleReal)param_.getValue("transformation_space:max_shift:MZ");
	
	      shift_bounding_box_.setMin(min);
	      shift_bounding_box_.setMax(max);
	
	      min[RawDataPoint2D::RT] = (DoubleReal)param_.getValue("transformation_space:min_scaling:RT");
	      min[RawDataPoint2D::MZ] = (DoubleReal)param_.getValue("transformation_space:min_scaling:MZ");
	      max[RawDataPoint2D::RT] = (DoubleReal)param_.getValue("transformation_space:max_scaling:RT");
	      max[RawDataPoint2D::MZ] = (DoubleReal)param_.getValue("transformation_space:max_scaling:MZ");
	
	      scaling_bounding_box_.setMin(min);
	      scaling_bounding_box_.setMax(max);
	
	    }
	
	    /// To speed up the calculation of the final transformation, we confine the number of
	    /// considered point pairs. We match a point p in the model map only onto those points p'
	    /// in the scene map that lie in a certain mz intervall.
	    /// If  (p_mz - mz_bucket_size_) <= p'_mz <= (p_mz mz_bucket_size_) then p and p' are partners.
	    void preprocess_()
	    {
	    	//Select to 2000 most abundant data points only //Do this only once when the mode map is set
	      PeakPointerArray model_map(model_map_->begin(),model_map_->end());
	      model_map.sortByIntensity(true);
				if (model_map.size()>2000) model_map.resize(2000); //TODO make this a parameter
				model_map.sortByNthPosition(RawDataPoint2D::MZ);
				
	      PeakPointerArray scene_map(scene_map_->begin(),scene_map_->end());
	      scene_map.sortByIntensity(true);
	      if (scene_map.size()>2000) scene_map.resize(2000); //TODO make this a parameter
				scene_map.sortByNthPosition(RawDataPoint2D::MZ);
				
	      // for each element (rt_m,mz_m) of the model map
	      // search for corresponding elements (rt_i,mz_i) in the scene map
	      // which lie in a predefined mz intervall (mz_i in [mz_m-eps,mz_m+eps))
	
	      // take only elements of the model map which have partners in the scene map
	      typename PeakPointerArray::const_iterator it_first = scene_map.begin();
	      typename PeakPointerArray::const_iterator it_last = it_first;
	      for (UInt i = 0; i < model_map.size(); ++i)
	      {
	        DoubleReal act_mz = model_map[i].getMZ();
	        DoubleReal min_mz = act_mz - mz_bucket_size_;
	        DoubleReal max_mz = act_mz + mz_bucket_size_;
	
	        std::vector< const PointType* > partners;
	        // search for the left end of the intervall //TODO use binary search
	        while ( it_first!=scene_map.end()  && it_first->getMZ()<min_mz)
	        {
	          ++it_first;
	        }
	
	        it_last = it_first;
	
	        // search for the right end of the intervall
	        while ((it_last != scene_map.end())  && (it_last->getMZ() < max_mz))
	        {
	          ++it_last;
	        }
	
	        for (typename PeakPointerArray::const_iterator it = it_first; it != it_last; ++it)
	        {
	        	partners.push_back(&(*it));
	        }
	        
	        if (partners.size() > 0)
	        {
	          model_map_red_.push_back(model_map[i]);
	          scene_map_partners_.push_back(partners);
	        }
	      }
	
	      // Compute shift_bucket_size_ and num_buckets.
	      DPosition<2> diagonal_shift = shift_bounding_box_.diagonal();
	      DPosition<2> diagonal_scaling = scaling_bounding_box_.diagonal();
	
        num_buckets_shift_ = (int)(ceil(diagonal_shift[0]/shift_bucket_size_));
        num_buckets_scaling_ = (int)(ceil(diagonal_scaling[0]/scaling_bucket_size_));
        shift_bucket_size_ = diagonal_shift[0] / (DoubleReal)num_buckets_shift_;
        scaling_bucket_size_ = diagonal_scaling[0] / (DoubleReal)num_buckets_scaling_;
	
	    } // preprocess_
	
	    /// Compute the transformations between each point pair in the model map and each point pair in the scene map
	    /// and hash the affine transformation.
	    void hashAffineTransformations_(DoubleReal total_int_model_map, DoubleReal total_int_scene_map )
	    {
	      // take each point pair in the model map
	      UInt n = model_map_red_.size();
	      for (UInt i = 0; i < n; ++i)
	      {
	        // take only the next 10 neighbours in m/z as partner in the model map
	        //        UInt k=((i+50)>= n) ? n : (i+50);
	        for (UInt j = i+1; j < n; ++j)
	        {
	          // avoid cross mappings (i,j) -> (k,l) (e.g. i_rt < j_rt and k_rt > l_rt)
	          // and point pairs with equal retention times (e.g. i_rt == j_rt)
	          DPosition<2> diff = model_map_red_[i].getPosition() - model_map_red_[j].getPosition();
	          // and compute the affine transformation to all corresponding points pair in the scene map
	          std::vector< const PointType* >& partners_i = scene_map_partners_[i];
	          std::vector< const PointType* >& partners_j  = scene_map_partners_[j];
	          UInt m = partners_i.size();
	          UInt p = partners_j.size();
	
	          for (UInt k = 0; k < m; ++k)
	          {
	            for (UInt l = 0; l < p; ++l)
	            {
	              DPosition<2> diff_2 = (partners_j[l]->getPosition() - partners_i[k]->getPosition());
	
	              // compute the transformation (i,j) -> (k,l)
	              DPosition<2> shift;
	              DPosition<2> scaling;
	              bool transformation_ok[2];
	              transformation_ok[0] = false;
	              transformation_ok[1] = false;
	
	              if ((fabs(diff[RawDataPoint2D::RT]) > 0.001) && (fabs(diff_2[RawDataPoint2D::RT]) > 0.001))
	              {
	                scaling[RawDataPoint2D::RT] = (model_map_red_[j].getRT() - model_map_red_[i].getRT()) /(partners_j[l]->getRT() - partners_i[k]->getRT());
	                shift[RawDataPoint2D::RT] =  model_map_red_[i].getRT() - partners_i[k]->getRT()*scaling[RawDataPoint2D::RT];
	                transformation_ok[RawDataPoint2D::RT] = true;
	              }
	
	
	              if ((fabs(diff[RawDataPoint2D::MZ]) > 0.001) && (fabs(diff_2[RawDataPoint2D::MZ]) > 0.001))
	              {
	                scaling[RawDataPoint2D::MZ] = (model_map_red_[j].getMZ() - model_map_red_[i].getMZ()) /(partners_j[l]->getMZ() - partners_i[k]->getMZ());
	                shift[RawDataPoint2D::MZ] =  model_map_red_[i].getMZ() - partners_i[k]->getMZ()*scaling[RawDataPoint2D::MZ];
	                transformation_ok[RawDataPoint2D::MZ] = true;
	              }
	
	              // compute the hash indices
	              int bucket_shift_index;
	              int bucket_scaling_index;
	              DoubleReal bucket_fraction_shift = 0.;
	              DoubleReal bucket_fraction_scaling = 0.;
	
	              DoubleReal int_i = model_map_red_[i].getIntensity()/total_int_model_map;
	              DoubleReal int_k = partners_i[k]->getIntensity()/total_int_scene_map;
	              DoubleReal int_j = model_map_red_[j].getIntensity()/total_int_model_map;
	              DoubleReal int_l = partners_j[l]->getIntensity()/total_int_scene_map;
	              DoubleReal int_similarity_jl = (int_j < int_l) ? int_j/int_l : int_l/int_j;
	              DoubleReal int_similarity_ik = (int_i < int_k) ? int_i/int_k : int_k/int_i;
	              DoubleReal int_similarity =  int_similarity_ik + int_similarity_jl;
	
	              // check if the transformation paramerters lie in the hash map
	              if ( transformation_ok[RawDataPoint2D::RT] && (scaling_bounding_box_.min()[RawDataPoint2D::RT] < scaling[RawDataPoint2D::RT]) && (scaling_bounding_box_.max()[RawDataPoint2D::RT] > scaling[RawDataPoint2D::RT])
	                   && (shift_bounding_box_.min()[RawDataPoint2D::RT] < shift[RawDataPoint2D::RT]) && (shift_bounding_box_.max()[RawDataPoint2D::RT] > shift[RawDataPoint2D::RT]))
	              {
	                bucket_fraction_shift = (shift[RawDataPoint2D::RT] - shift_bounding_box_.min()[RawDataPoint2D::RT]) / shift_bucket_size_;  // floating point division
	                bucket_fraction_scaling = (scaling[RawDataPoint2D::RT] - scaling_bounding_box_.min()[RawDataPoint2D::RT]) / scaling_bucket_size_;  // floating point division
	                bucket_shift_index    = (int) bucket_fraction_shift; // round down (yes we are >= 0)
	                bucket_scaling_index    = (int) bucket_fraction_scaling; // round down (yes we are >= 0)
	                bucket_fraction_shift -= bucket_shift_index;          // fractional part
	                bucket_fraction_scaling -= bucket_scaling_index;          // fractional part
	
	                // hash the transformation if possible
	                DoubleReal bucket_fraction_complement_shift = 1.;
	                DoubleReal bucket_fraction_complement_scaling = 1.;
	                bucket_fraction_complement_shift -= bucket_fraction_shift;
	                bucket_fraction_complement_scaling -= bucket_fraction_scaling;
	
	                // Distribute the vote of the shift among the four neighboring buckets.
	                DoubleReal factor;
	                // rt
	                factor = bucket_fraction_complement_shift * bucket_fraction_complement_scaling * int_similarity;
	                rt_hash_[PairType(bucket_shift_index, bucket_scaling_index)] += factor;
	
	                factor = bucket_fraction_complement_shift * bucket_fraction_scaling * int_similarity;
	                rt_hash_[PairType(bucket_shift_index, bucket_scaling_index + 1 )] += factor;
	
	                factor = bucket_fraction_shift * bucket_fraction_complement_scaling * int_similarity;
	                rt_hash_[PairType( bucket_shift_index + 1, bucket_scaling_index )] += factor;
	
	                factor = bucket_fraction_shift * bucket_fraction_scaling * int_similarity_ik;
	                rt_hash_[PairType( bucket_shift_index + 1, bucket_scaling_index + 1 )] += factor;
	              }
	            } // for l
	          } // for k
	        } // for j
	      } // for i
	
	    } // hashAffineTransformations_
	
	
	    /// After the hashing phase, the best transformation, that is the transformation with the most votes is determined.
	    void estimateFinalAffineTransformation_(LinearMapping& mapping)
	    {
	      // search for the maximal vote parameter of the rt transformation
	      PairType max_element_index_rt;
	      DoubleReal act_max_rt = 0;
	      for (typename AffineTransformationMapType::const_iterator it = rt_hash_.begin(); it != rt_hash_.end(); ++it)
	      {
	        if (it->second > act_max_rt)
	        {
	          max_element_index_rt = it->first;
	          act_max_rt = it->second;
	        }
	      }
	
	      // Compute a weighted average of the transformation parameters nearby the max_element_index.
	      DPosition<2> rt_trafo;
	      DPosition<2> rt_bounding_box_min(shift_bounding_box_.min()[RawDataPoint2D::RT],scaling_bounding_box_.min()[RawDataPoint2D::RT]);
	      int rt_run_indices[2];
	      DoubleReal quality = 0;
	      
	      for ( rt_run_indices[0]  = std::max ( int (max_element_index_rt.first - bucket_window_shift_), 0 );
	            rt_run_indices[0] <= std::min ( int (max_element_index_rt.first + bucket_window_shift_), num_buckets_shift_-1 );
	            ++rt_run_indices[0])
	      {
	        for ( rt_run_indices[1]  = std::max ( int (max_element_index_rt.second - bucket_window_scaling_), 0 );
	              rt_run_indices[1] <= std::min ( int (max_element_index_rt.second + bucket_window_scaling_), num_buckets_scaling_-1 );
	              ++rt_run_indices[1])
	
	        {
	          DPosition<2> contribution_position(shift_bucket_size_,scaling_bucket_size_);
	          // is the neighbouring bucket in the map?
	          typename AffineTransformationMapType::const_iterator it = rt_hash_.find(PairType(rt_run_indices[0], rt_run_indices[1]));
	          if ( it != rt_hash_.end())
	          {
	            for ( UInt dimension = 0; dimension < 2; ++dimension)
	            {
	              contribution_position[dimension] *= rt_run_indices[dimension];
	            }
	            contribution_position += rt_bounding_box_min;
	            DoubleReal contribution_quality = it->second;
	            quality += contribution_quality;
	            contribution_position *= contribution_quality;
	            rt_trafo += contribution_position;
	          }
	        }
	      }
	
	      if ( quality != 0 )
	      {
	        rt_trafo /= quality;
	      }
	
	      mapping.setSlope(rt_trafo[1]);
	      mapping.setIntercept(rt_trafo[0]);
	
			}
	
	    /// Reduced model map which contains only elements of the model map which have a partner in the scene map
	    ElementMapType model_map_red_;
	
	    /// Partner elements in the scene map
	    std::vector< std::vector< const PointType* > > scene_map_partners_;
	
	    /// Hash map of all transformations in the rt dimension (contains the affine transformation parameters)
	    AffineTransformationMapType rt_hash_;
	
	    /// Hash map of all transformations in the mz dimension (contains the affine transformation parameters)
	    AffineTransformationMapType mz_hash_;
	
	    /// Bounding box of the shift parameters
	    DBoundingBox<2> shift_bounding_box_;
	
	    /// Bounding box of the scaling parameters
	    DBoundingBox<2> scaling_bounding_box_;
	
	    /// Diagonal size of each shift bucket
	    DoubleReal shift_bucket_size_;
	
	    /// Diagonal size of each scaling bucket
	    DoubleReal scaling_bucket_size_;
	
	    /// Number of shift buckets
	    int num_buckets_shift_;
	
	    /// Number of scaling buckets
	    int num_buckets_scaling_;
	
	    /// Number of neighbouring shift buckets to be considered computing the final transformations
	    UInt bucket_window_shift_;
	
	    /// Number of neighbouring scaling buckets to be considered computing the final transformations
	    UInt bucket_window_scaling_;
	
	    /// Maximum deviation in mz of two partner points
	    DoubleReal mz_bucket_size_;
    
  }
  ; // PoseClusteringAffineSuperimposer
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_PoseClusteringAffineSuperimposer_H
