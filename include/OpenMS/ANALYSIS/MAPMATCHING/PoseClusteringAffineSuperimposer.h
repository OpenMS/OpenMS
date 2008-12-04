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

	It works on two element maps and computes an affine transformation that maps
	the elements of second map as near as possible to the elements in the first
	map.

	This superimposer hashes affine transformations between pairs of features in
	map one and pairs of features in map two.  Then, it finds the transformation
	with most votes.

	@todo Do all the todos in the code, then move stuff to .C file (Clemens)

	@htmlinclude OpenMS_PoseClusteringAffineSuperimposer.parameters

	@ingroup MapAlignment
	*/
	class PoseClusteringAffineSuperimposer
		: public BaseSuperimposer
	{
	 public:

		/// Constructor
		PoseClusteringAffineSuperimposer()
			: BaseSuperimposer()
		{
			setName(getProductName());

			// TODO rename mz_bucket_size to something more appropriate (mz_pair_max_distance?)
			defaults_.setValue("mz_bucket_size",0.5,"Maximum of m/z deviation of corresponding elements in different maps.");
			// defaults_.setValue("rt_pair_min_distance",10,"Minimum difference of retention time for pairs");
			defaults_.setValue("num_used_points",2000,"The maximum number of points used.  Points are selected by intensity.");
			defaults_.setValue("shift_bucket_size",10.0,"Defines the shift parameter's bucket size during histograming.",StringList::create("advanced"));
			defaults_.setValue("scaling_bucket_size",0.01,"Defines the scaling parameter's bucket size during histograming.",StringList::create("advanced"));
			defaults_.setValue("bucket_window_shift",2,"Number of surrounding buckets of element indices to be considered when computing the shift parameter.",StringList::create("advanced"));
			defaults_.setValue("bucket_window_scaling",2,"Number of surrounding buckets of element indices to be considered when computing the scaling parameter.",StringList::create("advanced"));
			defaults_.setValue("min_shift",-1000.0,"Minimal shift parameter which is considered during histogramming.",StringList::create("advanced"));
			defaults_.setValue("max_shift",1000.0,"Maximal shift parameter which is considered during histogramming.",StringList::create("advanced"));
			defaults_.setValue("min_scaling",0.5,"Minimal scaling parameter which is considered during histogramming.",StringList::create("advanced"));
			defaults_.setValue("max_scaling",2.0,"Maximal scaling parameter which is considered during histogramming.",StringList::create("advanced"));

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
		virtual void run ( const std::vector<ConsensusMap>& maps,
											 std::vector<TransformationDescription>& transformations
										 )
		{
			if (maps.size()!=2) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Excactly two input maps are required");

			// clear the hash table for parameters of affine transformation
			rt_hash_.clear();

			startProgress(0,120,"affine pose clustering");
			setProgress(0);
			//************************************************************************************
			// Preprocessing

			// To speed up the calculation of the final transformation, we confine the number of
			// considered point pairs. We match a point p in the model map only onto those points p'
			// in the scene map that lie in a certain mz intervall.
			// If  (p_mz - mz_bucket_size_) <= p'_mz <= (p_mz + mz_bucket_size_) then p and p' are partners.

			//Select the most abundant data points only
			UInt num_used_points = param_.getValue("num_used_points");

			PeakPointerArray model_map_ini(maps[0].begin(),maps[0].end());
			model_map_ini.sortByIntensity(true);
			if (model_map_ini.size()>num_used_points) model_map_ini.resize(num_used_points);
			model_map_ini.sortByComparator(Peak2D::MZLess());

			setProgress(2);

			PeakPointerArray scene_map_ini(maps[1].begin(),maps[1].end());
			scene_map_ini.sortByIntensity(true);
			if (scene_map_ini.size()>num_used_points) scene_map_ini.resize(num_used_points);
			scene_map_ini.sortByComparator(Peak2D::MZLess());

			setProgress(4);

			// No more changes after this point (we tend to have annoying issues
			// with const_iterator versus iterator)
			const PeakPointerArray & model_map(model_map_ini);
			const PeakPointerArray & scene_map(scene_map_ini);

			// compute total intensities of both maps for normalisation
			DoubleReal total_int_model_map = 0;
			for (UInt i = 0; i < model_map.size(); ++i)
			{
				total_int_model_map += model_map[i].getIntensity();
			}

			setProgress(6);

			DoubleReal total_int_scene_map = 0;
			for (UInt i = 0; i < scene_map.size(); ++i)
			{
				total_int_scene_map += scene_map[i].getIntensity();
			}

			setProgress(8);

			// Compute shift_bucket_size_ and num_buckets.
			DPosition<1> diagonal_shift = shift_bounding_box_.diagonal();
			DPosition<1> diagonal_scaling = scaling_bounding_box_.diagonal();

			num_buckets_shift_ = (int)(ceil(diagonal_shift[0]/shift_bucket_size_));
			num_buckets_scaling_ = (int)(ceil(diagonal_scaling[0]/scaling_bucket_size_));
			shift_bucket_size_ = diagonal_shift[0] / (DoubleReal)num_buckets_shift_;
			scaling_bucket_size_ = diagonal_scaling[0] / (DoubleReal)num_buckets_scaling_;

			setProgress(10);

			//************************************************************************************
			// Hashing

			// Compute the transformations between each point pair in the model map
			// and each point pair in the scene map and hash the affine
			// transformation.

			UInt const model_map_size = model_map.size(); // i j
			UInt const scene_map_size = scene_map.size(); // k l

			// first point in model map
			for ( UInt i = 0, k_low = 0, k_high = 0;
						i < model_map_size-1;
						++i
					)
			{
				setProgress(10+Real(i)/model_map_size*100.f);

				// Adjust window in scene map
				while ( k_low < scene_map_size &&
								scene_map[k_low].getMZ() < model_map[i].getMZ() - mz_bucket_size_
							) ++k_low ;
				while ( k_high < scene_map_size &&
								scene_map[k_high].getMZ() <= model_map[i].getMZ() + mz_bucket_size_
							) ++k_high ;

				// first point in scene map
				for ( UInt k = k_low; k < k_high; ++k )
				{
					// second point in model map
					for (UInt j = i+1, l_low = k_low, l_high = k_high; j < model_map_size; ++j)
					{
						// diff in model map
						DoubleReal diff_model = model_map[j].getRT() - model_map[i].getRT();
						// TODO magic alert!  use rt_pair_min_distance instead of 0.001
						if (fabs(diff_model) < 0.001) continue;

						// Adjust window in scene map
						while ( l_low < scene_map_size &&
										scene_map[l_low].getMZ() < model_map[j].getMZ() - mz_bucket_size_
									) ++l_low ;
						while ( l_high < scene_map_size &&
										scene_map[l_high].getMZ() <= model_map[j].getMZ() + mz_bucket_size_
									) ++l_high ;

						// second point in scene map
						for ( UInt l = l_low; l < l_high; ++l )
						{
							// diff in scene map
							DoubleReal diff_scene = scene_map[l].getRT() - scene_map[k].getRT();

							// avoid cross mappings (i,j) -> (k,l) (e.g. i_rt < j_rt and k_rt > l_rt)
							// and point pairs with equal retention times (e.g. i_rt == j_rt)
							// TODO magic alert!  use rt_pair_min_distance instead of 0.001
							if ( fabs(diff_scene) <= 0.001 || ((diff_model>0)!=(diff_scene>0)) ) continue;

							// compute the transformation (i,j) -> (k,l)
							DoubleReal scaling = diff_model / diff_scene;
							DoubleReal shift =  model_map[i].getRT() - scene_map[k].getRT()*scaling;

							DoubleReal int_i = model_map[i].getIntensity()/total_int_model_map;
							DoubleReal int_k = scene_map[k].getIntensity()/total_int_scene_map;
							DoubleReal int_j = model_map[j].getIntensity()/total_int_model_map;
							DoubleReal int_l = scene_map[l].getIntensity()/total_int_scene_map;
							DoubleReal int_similarity_jl = (int_j < int_l) ? int_j/int_l : int_l/int_j;
							DoubleReal int_similarity_ik = (int_i < int_k) ? int_i/int_k : int_k/int_i;
							DoubleReal int_similarity =  int_similarity_ik + int_similarity_jl;

							// check if the transformation parameters lie in the hash map
							if (
									( scaling_bounding_box_.min()[0] >= scaling ) ||
									( scaling_bounding_box_.max()[0] <= scaling ) ||
									( shift_bounding_box_.min()[0] >= shift ) ||
									( shift_bounding_box_.max()[0] <= shift )
								 ) continue;

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

						} // l
					} // j
				} // k
			} // i 

			setProgress(110);

			//************************************************************************************
			// Estimate transform

			// search for the maximal vote parameter
			PairType max_element_index_rt;
			DoubleReal act_max_rt = 0;
			for (AffineTransformationMapType::const_iterator it = rt_hash_.begin(); it != rt_hash_.end(); ++it)
			{
				if (it->second > act_max_rt)
				{
					max_element_index_rt = it->first;
					act_max_rt = it->second;
				}
			}

			setProgress(112);

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
					AffineTransformationMapType::const_iterator it = rt_hash_.find(PairType(shift_index, scale_index));
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

			setProgress(116);

			// set trafo
			transformations.clear();
			transformations.resize(1);
			transformations[0].setName("linear");
			transformations[0].setParam("intercept",rt_shift / quality_sum * shift_bucket_size_ + shift_bounding_box_.min()[0]);
			transformations[0].setParam("slope",rt_scale / quality_sum * scaling_bucket_size_ + scaling_bounding_box_.min()[0]);

			setProgress(120);
			endProgress();

			return;
		}

		/// Returns an instance of this class
		static BaseSuperimposer* create()
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
		typedef ConstRefVector<ConsensusMap> PeakPointerArray;
		typedef std::pair<int,int> PairType;
		typedef std::map< PairType, DoubleReal> AffineTransformationMapType;
		//@}

		//docu in base class
		virtual void updateMembers_()
		{
			mz_bucket_size_ = param_.getValue("mz_bucket_size");
			shift_bucket_size_ = param_.getValue("shift_bucket_size");
			scaling_bucket_size_ = param_.getValue("scaling_bucket_size");
			bucket_window_shift_ = param_.getValue("bucket_window_shift");
			bucket_window_scaling_ = param_.getValue("bucket_window_scaling");

			DPosition<1> min;
			DPosition<1> max;
			min[0] = param_.getValue("min_shift");
			max[0] = param_.getValue("max_shift");

			shift_bounding_box_.setMin(min);
			shift_bounding_box_.setMax(max);

			min[0] = param_.getValue("min_scaling");
			max[0] = param_.getValue("max_scaling");

			scaling_bounding_box_.setMin(min);
			scaling_bounding_box_.setMax(max);

		}

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

#endif  // OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGAFFINESUPERIMPOSER_H
