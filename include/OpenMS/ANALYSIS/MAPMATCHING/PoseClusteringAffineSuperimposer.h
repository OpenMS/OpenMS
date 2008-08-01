// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H

#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/ConstRefVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>

#include <fstream>
#include <vector>
#include <map>
#include <math.h>

namespace OpenMS
{

  /**
    @brief Superimposer that uses a voting scheme to find a good affine transformation.

	  It works on two element maps and computes an affine transformation, that maps the 
	  elements of second map as near as possible to the elements in the first map.
	  An element can be a DPeak, a DFeature or ConsensusFeature.
	
	  This superimposer hashes affine transformations between pairs of features in map one and pairs of features in map two.
	  Then, it finds the transformation with most votes.
		
		@todo Do all the todos in the code (Clemens)
		
		@ref PoseClusteringAffineSuperimposer_Parameters are explained on a separate page.        

		@ingroup MapAlignment
  */
  template < typename MapT = FeatureMap<> >
  class PoseClusteringAffineSuperimposer 
  	: public BaseSuperimposer< MapT >
  {
  	public:
  		
  		///Base class type definition
			typedef BaseSuperimposer< MapT > Base;
			///Input map type
	    typedef typename Base::ElementMapType ElementMapType;
	    
	    using Base::param_;
	    using Base::defaultsToParam_;
	    using Base::defaults_;

	    /// Constructor
	    PoseClusteringAffineSuperimposer()
	    	: Base()
	    {
	      Base::setName(getProductName());
	
	      defaults_.setValue("mz_bucket_size",0.5,"An estimate of m/z deviation of corresponding elements in different maps.");
	      defaults_.setValue("num_used_points",2000,"The number of points used.\nThe most intense points are used");
	      defaults_.setValue("shift_bucket_size",10.0,"Defines the shift parameter's bucket size during histograming.",true);
	      defaults_.setValue("scaling_bucket_size",0.01,"Defines the scaling parameter's bucket size during histograming.",true);
	      defaults_.setValue("bucket_window_shift",2,"Number of surrounding buckets of element indices to be considered when computing the shift parameter.",true);
	      defaults_.setValue("bucket_window_scaling",2,"Number of surrounding buckets of element indices to be considered when computing the scaling parameter.",true);
	      defaults_.setValue("min_shift",-1000.0,"Minimal shift parameter which is considered during histogramming.",true);
	      defaults_.setValue("max_shift",1000.0,"Maximal shift parameter which is considered during histogramming.",true);
	      defaults_.setValue("min_scaling",0.5,"Minimal scaling parameter which is considered during histogramming.",true);
	      defaults_.setValue("max_scaling",2.0,"Maximal scaling parameter which is considered during histogramming.",true);
	
	      defaultsToParam_();
	    }

	    /// Destructor
	    virtual ~PoseClusteringAffineSuperimposer()
	    {
	    }

	    /**
	    	@brief Estimates the transformation and fills the given mapping function
	    	
	    	@note Exactly two input maps must be given.
	    	
	    	@exception IllegalArgument is thrown if the input maps are invalid.
	    */
	    virtual void run(const std::vector<ElementMapType>& maps, std::vector<TransformationDescription>& transformations)
	    {
	    	if (maps.size()!=2) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Excactly two input maps are required");
	    	
        // clear 
        scene_map_partners_.clear();
        model_map_red_.clear();
        rt_hash_.clear();

				//************************************************************************************
				// Preprocessing
	
		    // To speed up the calculation of the final transformation, we confine the number of
		    // considered point pairs. We match a point p in the model map only onto those points p'
		    // in the scene map that lie in a certain mz intervall.
		    // If  (p_mz - mz_bucket_size_) <= p'_mz <= (p_mz mz_bucket_size_) then p and p' are partners.

	    	//Select to 2000 most abundant data points only //Do this only once when the mode map is set
	      UInt num_used_points = param_.getValue("num_used_points");
	      
	      PeakPointerArray model_map(maps[0].begin(),maps[0].end());
	      model_map.sortByIntensity(true);
				if (model_map.size()>num_used_points) model_map.resize(num_used_points);
				model_map.sortByNthPosition(Peak2D::MZ);
				
	      PeakPointerArray scene_map(maps[1].begin(),maps[1].end());
	      scene_map.sortByIntensity(true);
	      if (scene_map.size()>num_used_points) scene_map.resize(num_used_points);
				scene_map.sortByNthPosition(Peak2D::MZ);

        // compute total intensities of both maps for normalisation
        DoubleReal total_int_model_map = 0;
        for (UInt i = 0; i < model_map.size(); ++i)
        {
          total_int_model_map += model_map[i].getIntensity();
        }

        DoubleReal total_int_scene_map = 0;
        for (UInt i = 0; i < scene_map.size(); ++i)
        {
          total_int_scene_map += scene_map[i].getIntensity();
        }

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
	        	partners.push_back(&(*it)); //TODO do not store the parterns, but do what has to be done right here
	        }
	        
	        if (partners.size() > 0)
	        {
	          model_map_red_.push_back(model_map[i]);
	          scene_map_partners_.push_back(partners);
	        }
	      }
	
	      // Compute shift_bucket_size_ and num_buckets.
	      DPosition<1> diagonal_shift = shift_bounding_box_.diagonal();
	      DPosition<1> diagonal_scaling = scaling_bounding_box_.diagonal();
	
        num_buckets_shift_ = (int)(ceil(diagonal_shift[0]/shift_bucket_size_));
        num_buckets_scaling_ = (int)(ceil(diagonal_scaling[0]/scaling_bucket_size_));
        shift_bucket_size_ = diagonal_shift[0] / (DoubleReal)num_buckets_shift_;
        scaling_bucket_size_ = diagonal_scaling[0] / (DoubleReal)num_buckets_scaling_;

				//************************************************************************************
				// Hashing
	    	
	    	// Compute the transformations between each point pair in the model map and each point pair in the scene map
	    	// and hash the affine transformation.
	      
	      // take each point pair in the model map
	      for (UInt i = 0; i < model_map_red_.size(); ++i)
	      {
	        // take only the next 10 neighbours in m/z as partner in the model map
	        //        UInt k=((i+50)>= n) ? n : (i+50);
	        for (UInt j = i+1; j < model_map_red_.size(); ++j)
	        {
	        	//diff in model map
	          DoubleReal diff_model = model_map_red_[j].getRT() - model_map_red_[i].getRT();
	          if (fabs(diff_model) > 0.001)
	          {
		          // and compute the affine transformation to all corresponding points pair in the scene map
		          std::vector< const PointType* >& partners_i = scene_map_partners_[i];
		          std::vector< const PointType* >& partners_j  = scene_map_partners_[j];
		
		          for (UInt k = 0; k < partners_i.size(); ++k)
		          {
		            for (UInt l = 0; l < partners_j.size(); ++l)
		            {
		            	//diff in scene map
		              DoubleReal diff_scene = partners_j[l]->getRT() - partners_i[k]->getRT();
			          	
			          	// avoid cross mappings (i,j) -> (k,l) (e.g. i_rt < j_rt and k_rt > l_rt)
	          			// and point pairs with equal retention times (e.g. i_rt == j_rt)
		              if (fabs(diff_scene) > 0.001 && ((diff_model>0)==(diff_scene>0)))
		              {
		              	// compute the transformation (i,j) -> (k,l)
		                DoubleReal scaling = diff_model / diff_scene;
		                DoubleReal shift =  model_map_red_[i].getRT() - partners_i[k]->getRT()*scaling;
		
			              DoubleReal int_i = model_map_red_[i].getIntensity()/total_int_model_map;
			              DoubleReal int_k = partners_i[k]->getIntensity()/total_int_scene_map;
			              DoubleReal int_j = model_map_red_[j].getIntensity()/total_int_model_map;
			              DoubleReal int_l = partners_j[l]->getIntensity()/total_int_scene_map;
			              DoubleReal int_similarity_jl = (int_j < int_l) ? int_j/int_l : int_l/int_j;
			              DoubleReal int_similarity_ik = (int_i < int_k) ? int_i/int_k : int_k/int_i;
			              DoubleReal int_similarity =  int_similarity_ik + int_similarity_jl;
			
			              // check if the transformation paramerters lie in the hash map
			              if ( (scaling_bounding_box_.min()[0] < scaling) && (scaling_bounding_box_.max()[0] > scaling)
			                   && (shift_bounding_box_.min()[0] < shift) && (shift_bounding_box_.max()[0] > shift))
			              {
			                DoubleReal bucket_fraction_shift = (shift - shift_bounding_box_.min()[0]) / shift_bucket_size_;  // floating point division
			                DoubleReal bucket_fraction_scaling = (scaling - scaling_bounding_box_.min()[0]) / scaling_bucket_size_;  // floating point division
			                Int bucket_shift_index    = (int) bucket_fraction_shift; // round down (yes we are >= 0)
			                Int bucket_scaling_index    = (int) bucket_fraction_scaling; // round down (yes we are >= 0)
			                bucket_fraction_shift -= bucket_shift_index;          // fractional part
			                bucket_fraction_scaling -= bucket_scaling_index;          // fractional part
			
			                // hash the transformation if possible
			                DoubleReal bucket_fraction_complement_shift = 1.;
			                DoubleReal bucket_fraction_complement_scaling = 1.;
			                bucket_fraction_complement_shift -= bucket_fraction_shift;
			                bucket_fraction_complement_scaling -= bucket_fraction_scaling;
			
			                // Distribute the vote of the shift among the four neighboring buckets.
			                rt_hash_[PairType(bucket_shift_index, bucket_scaling_index)] += bucket_fraction_complement_shift * bucket_fraction_complement_scaling * int_similarity;
			                rt_hash_[PairType(bucket_shift_index, bucket_scaling_index + 1 )] += bucket_fraction_complement_shift * bucket_fraction_scaling * int_similarity;
			                rt_hash_[PairType( bucket_shift_index + 1, bucket_scaling_index )] += bucket_fraction_shift * bucket_fraction_complement_scaling * int_similarity;
			                rt_hash_[PairType( bucket_shift_index + 1, bucket_scaling_index + 1 )] += bucket_fraction_shift * bucket_fraction_scaling * int_similarity;
			              }
		              }
		            }
	            }
	          }
	        }
	      }

				//************************************************************************************
				// Estimate transform
	      
	      // search for the maximal vote parameter
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
	      DoubleReal rt_shift = 0.0;
	      DoubleReal rt_scale = 0.0;
	      Int shift_index=0;
	      Int scale_index=0;
	      DoubleReal quality_sum = 0.0;
	      for ( shift_index  = std::max ( int (max_element_index_rt.first - bucket_window_shift_), 0 );
	            shift_index <= std::min ( int (max_element_index_rt.first + bucket_window_shift_), num_buckets_shift_-1 );
	            ++shift_index)
	      {
	        for ( scale_index  = std::max ( int (max_element_index_rt.second - bucket_window_scaling_), 0 );
	              scale_index <= std::min ( int (max_element_index_rt.second + bucket_window_scaling_), num_buckets_scaling_-1 );
	              ++scale_index)
	
	        {
	          // TODO use lower bound instead of find
	          typename AffineTransformationMapType::const_iterator it = rt_hash_.find(PairType(shift_index, scale_index));
	          if ( it != rt_hash_.end())
	          {
	            quality_sum += it->second;
	            rt_shift += shift_index * it->second;
	            rt_scale += scale_index * it->second;
	          }
	        }
	      }
	
	      if ( quality_sum == 0.0 )
	      {
	      	throw Exception::DivisionByZero(__FILE__,__LINE__,__PRETTY_FUNCTION__);
	      }
				
				//set trafo
				transformations.clear();
				transformations.resize(1);
				transformations[0].setName("linear");
	      transformations[0].setParam("intercept",rt_shift / quality_sum * shift_bucket_size_ + shift_bounding_box_.min()[0]);
	      transformations[0].setParam("slope",rt_scale / quality_sum * scaling_bucket_size_ + scaling_bounding_box_.min()[0]);
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
			///@name Internal type definitions
			//@{
	    typedef typename ElementMapType::value_type PointType;
	    typedef ConstRefVector<ElementMapType> PeakPointerArray;
	    typedef std::pair<int,int> PairType;
	    typedef std::map< PairType, DoubleReal> AffineTransformationMapType;
			//@}

	  	//docu in base class
	    virtual void updateMembers_()
	    {
	      mz_bucket_size_ = (DoubleReal)param_.getValue("mz_bucket_size");
	      shift_bucket_size_ = (DoubleReal)param_.getValue("shift_bucket_size");
	      scaling_bucket_size_ = (DoubleReal)param_.getValue("scaling_bucket_size");
	      bucket_window_shift_  = (UInt)param_.getValue("bucket_window_shift");
	      bucket_window_scaling_ = (UInt)param_.getValue("bucket_window_scaling");
	
	      DPosition<1> min;
	      DPosition<1> max;
	      min[0] = (DoubleReal)param_.getValue("min_shift");
	      max[0] = (DoubleReal)param_.getValue("max_shift");
	
	      shift_bounding_box_.setMin(min);
	      shift_bounding_box_.setMax(max);
	
	      min[0] = (DoubleReal)param_.getValue("min_scaling");
	      max[0] = (DoubleReal)param_.getValue("max_scaling");
	
	      scaling_bounding_box_.setMin(min);
	      scaling_bounding_box_.setMax(max);
	
	    }
	
	    /// Reduced model map which contains only elements of the model map which have a partner in the scene map
	    ElementMapType model_map_red_;
	
	    /// Partner elements in the scene map
	    std::vector< std::vector< const PointType* > > scene_map_partners_;
	
	    /// Hash map of all transformations in the rt dimension (contains the affine transformation parameters)
	    AffineTransformationMapType rt_hash_;
	
	    /// Bounding box of the shift parameters
	    DBoundingBox<1> shift_bounding_box_;
	
	    /// Bounding box of the scaling parameters
	    DBoundingBox<1> scaling_bounding_box_;
	
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
    
  };
  
} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_PoseClusteringAffineSuperimposer_H
