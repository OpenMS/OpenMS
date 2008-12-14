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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGSHIFTSUPERIMPOSER_H
#define OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGSHIFTSUPERIMPOSER_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/DATASTRUCTURES/Matrix.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>


#include <fstream>


namespace OpenMS
{
  /**
		@brief Superimposer that uses a voting scheme to find a good translation.
		
		It works similarly to the @ref PoseClusteringAffineSuperimposer , but only shifts are
		considered. Scaling is not corrected.
		
		@todo remove m/z correction and parameters; consider full rewrite (Clemens)
		
		@htmlinclude OpenMS_PoseClusteringShiftSuperimposer.parameters    

		@ingroup MapAlignment
  */
  class OPENMS_DLLAPI PoseClusteringShiftSuperimposer
		: public BaseSuperimposer
  {
  	public:

	    /// Constructor
	    PoseClusteringShiftSuperimposer()
	    	: BaseSuperimposer()
	    {
				setName(getProductName());
				
	      defaults_.setValue("input_map:bucket_size:RT",150.0,"Number of surrounding buckets of element indices to be considered when computing shifts.",StringList::create("advanced"));
	      defaults_.setValue("input_map:bucket_size:MZ",4.0,"Number of surrounding buckets of element indices to be considered when computing shifts.",StringList::create("advanced"));
	      defaults_.setValue("transformation_space:shift_bucket_size:RT",5.0,"Defines the shift parameter's bucket size during histograming.");
	      defaults_.setValue("transformation_space:shift_bucket_size:MZ",0.1,"Defines the shift parameter's bucket size during histograming.");
	      defaults_.setValue("input_map:bucket_window:RT",2,"Number of surrounding buckets of element indices to be considered when computing shifts.",StringList::create("advanced"));
	      defaults_.setValue("input_map:bucket_window:MZ",1,"Number of surrounding buckets of element indices to be considered when computing shifts.",StringList::create("advanced"));
	      defaults_.setValue("transformation_space:bucket_window_shift:RT",2,"Number of surrounding buckets of shift indices to be considered when computing shifts.",StringList::create("advanced"));
	      defaults_.setValue("transformation_space:bucket_window_shift:MZ",1,"Number of surrounding buckets of shift indices to be considered when computing shifts.",StringList::create("advanced"));
				subsections_.push_back("debug");
				
	      defaultsToParam_();
	    }
	
	    /// Destructor
	    virtual ~PoseClusteringShiftSuperimposer()
	    {
	    }
	
	    /**
	    	@brief Estimates the transformation and fills the given mapping function
	    	
	    	@note Exactly two input maps must be given.
	    	
	    	@exception IllegalArgument is thrown if the input maps are invalid.
	    */
	    virtual void run(const std::vector<ConsensusMap>& maps, std::vector<TransformationDescription>& transformations)
	    {
	    	if (maps.size()!=2) throw Exception::IllegalArgument(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Excactly two input maps are required");
	    	
	    	model_map_ = &(maps[MODEL]);
	    	scene_map_ = &(maps[SCENE]);
	    	
	      // clear the member
	      element_bucket_[Peak2D::RT].clear();
	      element_bucket_[Peak2D::MZ].clear();
	      shift_bucket_.clear();

				//clear trafo
				transformations.clear();
				transformations.resize(1);

	      computeElementBuckets_();
	      computeShiftBuckets_();
	      computeShift_(transformations[0]);
	    }
	
	    /// Returns an instance of this class
	    static BaseSuperimposer* create()
	    {
	      return new PoseClusteringShiftSuperimposer();
	    }
	
	    /// Returns the name of this module
	    static const String getProductName()
	    {
	      return "poseclustering_shift";
	    }
	    
	  protected:

	    /// Internal representation of a shift used in PoseClusteringShiftSuperimposer
	    struct Shift
	    {
				///Constructor
				Shift()
					: position(),
				    quality(0.0)
				{
				}
				
	      DPosition<2> position;
	      DoubleReal quality;
	    };

			///@name Type definitions
			//@{
	    typedef ConsensusMap::value_type PointType;
	    typedef Matrix < std::vector<UInt> > ElementBucketMatrixType;
	    typedef Matrix < DoubleReal > ShiftQualityMatrixType;
	    typedef Matrix < Shift > ShiftMatrixType;
			//@}
	  	
	  	///map index names
	  	enum
			{
				MODEL=0,
				SCENE=1
			};
	  	
	    virtual void updateMembers_()
	    {
	      shift_bucket_size_[0] = (DoubleReal)param_.getValue("transformation_space:shift_bucket_size:RT");
	      shift_bucket_size_[1] = (DoubleReal)param_.getValue("transformation_space:shift_bucket_size:MZ");
	      element_bucket_window_[0] = (UInt)param_.getValue("input_map:bucket_window:RT");
	      element_bucket_window_[1] = (UInt)param_.getValue("input_map:bucket_window:MZ");
	      shift_bucket_window_[0] = (UInt)param_.getValue("transformation_space:bucket_window_shift:RT");
	      shift_bucket_window_[1] = (UInt)param_.getValue("transformation_space:bucket_window_shift:MZ");
	      element_bucket_size_[0] = (DoubleReal)param_.getValue("input_map:bucket_size:RT");
	      element_bucket_size_[1] = (DoubleReal)param_.getValue("input_map:bucket_size:MZ");
	    }

	    /// Fill the buckets with the indices of the corresponding elements.
	    void computeElementBuckets_()
	    {
	      // Shorthands ...
	      DPosition<2> & fbs = element_bucket_size_;
	
				const ConsensusMap* map_array[2] = {model_map_, scene_map_};
	
	      for ( UInt map_index = 0; map_index < 2; ++map_index )
	      {
	        // Shorthands ...
	      	ConsensusMap const & fm = *(map_array[map_index]);
	        DBoundingBox<2> & fmpbb  = element_map_position_bounding_box_[map_index] ;
	        DBoundingBox<1> & fmibb  = element_map_intensity_bounding_box_[map_index];
	
	        fmpbb.clear();
	        fmibb.clear();
	
	        // Compute the bounding box for the element map, with respect to
	        // position and intensity.
	        for ( ConsensusMap::const_iterator fm_iter = fm.begin();
	              fm_iter != fm.end();
	              ++fm_iter
	            )
	        {
	          fmpbb.enlarge(fm_iter->getPosition());
	          fmibb.enlarge(fm_iter->getIntensity());
	        }
	      }
	
	      // Next we will enlarge each element_map_position_bounding_box_ such
	      // that all buckets will have the same diagonal.  To provide against
	      // rounding errors, we allocate one bucket more than needed (in each
	      // dimension) and shift the grid by one-half of the difference.
	      for ( UInt map_index = 0; map_index < 2; ++map_index )
	      {
	        // Shorthands ...
	        ConsensusMap const & fm = *(map_array[map_index]);
	        DBoundingBox<2> const & fmpbb  = element_map_position_bounding_box_[map_index] ;
	        DBoundingBox<2> & fmpbbe = element_map_position_bounding_box_enlarged_[map_index] ;
	        ElementBucketMatrixType       & fb     = element_bucket_[map_index];
	
	        // Compute num_buckets.  Compute extra margin to make bounding box a
	        // multiple of element buckets.
	        DPosition<2> const diagonal = fmpbb.diagonal();
	        DPosition<2> diagonal_enlarged;
	        int num_buckets[2];
	        for ( UInt dimension = 0; dimension < 2; ++dimension)
	        {
	          num_buckets[dimension] = int(1.1 + diagonal[dimension]/fbs[dimension]);
	          diagonal_enlarged[dimension] = fbs[dimension] * num_buckets[dimension];
	        }
	
	        // The extra margin.
	        DPosition<2> extra_element_bucket_size_(diagonal_enlarged-diagonal);
	        extra_element_bucket_size_ /= 2;
	
	        // Compute the enlarged element map bounding box accordingly.
	        fmpbbe.clear();
	        fmpbbe.enlarge( fmpbb.min() - extra_element_bucket_size_ );
	        fmpbbe.enlarge( fmpbb.max() + extra_element_bucket_size_ );
	
	        // Resize element_bucket_[map_index] accordingly.
	        fb.resize(num_buckets[Peak2D::RT],num_buckets[Peak2D::MZ]);
	
	        // Now, finally, we store the indices of the elements in their
	        // corresponding buckets.
	        DPosition<2> const & fmpbbe_min = fmpbbe.min();
	        for ( UInt index= 0; index < fm.size(); ++index )
	        {
	          DPosition<2> position = fm[index].getPosition() - fmpbbe_min;
	          fb ( UInt(position[Peak2D::RT]/fbs[Peak2D::RT]), UInt(position[Peak2D::MZ]/fbs[Peak2D::MZ]) ).push_back(index);
	        }
	
	        // Optionally, write debug output as specified in param.
	        if ( getParameters().exists("debug:feature_buckets_file") )
	        {
	          String const element_buckets_file = (String)getParameters().getValue("debug:feature_buckets_file") + String(map_index?"_SCENE":"_MODEL");
	          std::ofstream dump_file(element_buckets_file.c_str());
	          std::cerr << "### Writing "<<element_buckets_file<<std::endl;
	          dump_file << "# " << element_buckets_file << " generated " << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << std::endl;
	          dump_file << "# Positions of elements in non-empty element buckets" << std::endl;
	          for ( ElementBucketMatrixType::ConstIterator iter = fb.begin(); iter != fb.end(); ++iter)
	          {
	            if (iter->empty())
	              continue;
	            std::pair<UInt,UInt> row_col = fb.indexPair(iter-fb.begin());
	            dump_file << row_col.first << ' ' << row_col.second << " #bucket" << std::endl;
	            for ( std::vector<UInt>::const_iterator viter = iter->begin(); viter != iter->end(); ++viter)
	            {
	              dump_file << fm[*viter].getRT() <<' '<<fm[*viter].getMZ() << std::endl;
	            }
	            dump_file << std::endl;
	          }
	          dump_file << "# " << element_buckets_file << " EOF " << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << std::endl;
	        }
	      }
	    }
	
	    /**
	    	@brief Fill the buckets of shifts.
	
		    Note that computeElementBuckets_() must have been called before to make this work properly.
	    */
	    void computeShiftBuckets_()
	    {
	      // Shorthands ...
	      ShiftQualityMatrixType & tb = shift_bucket_;
	      DPosition<2> & tbs = shift_bucket_size_;
	      DBoundingBox<2> & tbb = shift_bounding_box_ ;
	      DBoundingBox<2> & tbbe = shift_bounding_box_enlarged_ ;
	      UInt const (&fbw)[2] = element_bucket_window_;
	
	      // Compute the bounding box for the shift map
	      {
	        tbb.clear();
	        tbb.enlarge ( element_map_position_bounding_box_[SCENE].min() - element_map_position_bounding_box_[MODEL].min() );
	        tbb.enlarge ( element_map_position_bounding_box_[SCENE].min() - element_map_position_bounding_box_[MODEL].max() );
	        tbb.enlarge ( element_map_position_bounding_box_[SCENE].max() - element_map_position_bounding_box_[MODEL].min() );
	        tbb.enlarge ( element_map_position_bounding_box_[SCENE].max() - element_map_position_bounding_box_[MODEL].max() );
	      }
	      // Next we will enlarge each bucket_size_ such that all buckets will
	      // have the same diagonal.  To provide against rounding errors, we
	      // allocate one bucket more than needed (in each dimension) and shift
	      // the grid by one-half.
	
	      DPosition<2> half_of_shift_bucket_size_(tbs);
	      half_of_shift_bucket_size_ /= 2;
	
	      // Adjust the enlarged shift map bounding box accordingly.
	      {
	        tbbe.clear();
	        tbbe.enlarge( tbb.min() - half_of_shift_bucket_size_ );
	        tbbe.enlarge( tbb.max() + half_of_shift_bucket_size_ );
	      }
	
	      // Compute shift_bucket_size_ and num_buckets.
	      DPosition<2> diagonal = tbbe.diagonal();
	      int num_buckets[2];
	      for ( UInt dimension = 0; dimension < 2; ++dimension)
	      {
	        num_buckets[dimension] = int(diagonal[dimension]/tbs[dimension]);
	        tbs[dimension] = diagonal[dimension] / num_buckets[dimension];
	      }
	
	      // Resize shift_bucket_ accordingly.
	      tb.resize(num_buckets[Peak2D::RT]+1,num_buckets[Peak2D::MZ]+1);
	
	      // Clear the shift buckets.
	      std::fill(tb.begin(),tb.end(),DoubleReal(0));
	
	
	      // Resize shift_matrix_ according to element_bucket_[Peak2D::MZ]
	      //         tm.resize(element_bucket_[Peak2D::MZ].sizePair());
	
	      // Now we store the shifts for all relevant element pairs in their
	      // corresponding buckets.  Each shift is distributed among its
	      // four neighboring "buckets", with weights according to the distances
	      // from these corner points.  Note that the outer two loops (over i and
	      // j) enumerate the "image" (element_bucket_[Peak2D::MZ]), then we search for
	      // "pre-images" (element_bucket_[0}) in the two inner loops (over k and
	      // l).  (And of course, finally, we enumerate all element pairs.)  This
	      // way we can associate the shifts vectors to buckets of the
	      // image, and when we will later apply it, we will not change the
	      // pre-image, which might be a consensus or so.
	
	      // progress dots
	      Int progress_dots = 0;
				if (this->param_.exists("debug::progress_dots"))
				{
	      	progress_dots = (Int)this->param_.getValue("debug:progress_dots");
			 	}
	
	      DPosition<2> const & tbbe_min = tbbe.min();
	
	      // Compute the index shift of corresponding element buckets of model and scene.
	      DPosition<2> const fmpbbe_min_offset =
	        element_map_position_bounding_box_enlarged_[SCENE].min() -
	        element_map_position_bounding_box_enlarged_[MODEL].min();
	      int const element_buckets_index_offset_RT = int ( fmpbbe_min_offset[Peak2D::RT] / element_bucket_size_[Peak2D::RT] );
	      int const element_buckets_index_offset_MZ = int ( fmpbbe_min_offset[Peak2D::MZ] / element_bucket_size_[Peak2D::MZ] );
	
	      // iterate over buckets of scene
	      for ( UInt scene_bucket_index_RT = 0;
	            scene_bucket_index_RT < element_bucket_[SCENE].rows();
	            ++scene_bucket_index_RT
	          )
	      {
	        for ( UInt scene_bucket_index_MZ = 0;
	              scene_bucket_index_MZ < element_bucket_[SCENE].cols();
	              ++scene_bucket_index_MZ
	            )
	        {
	
	          // compute the corresponding bucket in the model
	          int const model_bucket_index_center_RT = scene_bucket_index_RT + element_buckets_index_offset_RT;
	          int const model_bucket_index_center_MZ = scene_bucket_index_MZ + element_buckets_index_offset_MZ;
	
	          // iterate over buckets of model
	          for ( int model_bucket_index_RT
	                =  std::max<int>( model_bucket_index_center_RT - fbw[Peak2D::RT], 0 );
	                model_bucket_index_RT
	                <= std::min<int>( model_bucket_index_center_RT + fbw[Peak2D::RT], element_bucket_[MODEL].rows()-1 );
	                ++model_bucket_index_RT
	              )
	          {
	            for ( int model_bucket_index_MZ
	                  =  std::max<int>( model_bucket_index_center_MZ - fbw[Peak2D::MZ], 0 );
	                  model_bucket_index_MZ
	                  <= std::min<int>( model_bucket_index_center_MZ + fbw[Peak2D::MZ], element_bucket_[MODEL].cols()-1 );
	                  ++model_bucket_index_MZ
	                )
	            {
	              // iterate over pairs of elements for this pair of buckets
	              int number_of_considered_element_pairs_for_this_pair_of_buckets = 0;
	              std::vector<UInt> const & model_element_bucket
	              = element_bucket_[MODEL]
	                ( model_bucket_index_RT, model_bucket_index_MZ );
	              for ( std::vector<UInt>::const_iterator model_iter = model_element_bucket.begin();
	                    model_iter != model_element_bucket.end();
	                    ++model_iter
	                  )
	              {
	                std::vector<UInt> const & scene_element_bucket
	                = element_bucket_[SCENE]( scene_bucket_index_RT, scene_bucket_index_MZ );
	                for ( std::vector<UInt>::const_iterator scene_iter = scene_element_bucket.begin();
	                      scene_iter != scene_element_bucket.end();
	                      ++scene_iter
	                    )
	                {
	                  // Compute the shift corresponding to a pair of elements.
	                  Shift shift = shift_( (*model_map_)[*model_iter], (*scene_map_)[*scene_iter] );
	
	                  DPosition<2> tpwm = shift.position;
	                  tpwm -= tbbe_min;
	
	                  DoubleReal  const & tq = shift.quality;
	
	                  // Compute the bucket index (the lowest of the four) for
	                  // this shift.  Also compute the fractional part of
	                  // the position within the bucket.
	                  UInt bucket_index[2];
	                  DPosition<2> bucket_fraction;
	                  for ( UInt dimension = 0; dimension < 2; ++dimension )
	                  {
	                    bucket_fraction[dimension] = tpwm[dimension] / tbs[dimension];  // floating point division
	                    bucket_index[dimension]    = (UInt) bucket_fraction[dimension]; // round down (yes we are >= 0)
	                    bucket_fraction[dimension] -= bucket_index[dimension];          // fractional part
	                  }
	                  DPosition<2> bucket_fraction_complement(1,1);
	                  bucket_fraction_complement -= bucket_fraction;
	
	                  // Distribute the quality of the shift among the four neighboring buckets.
	                  DoubleReal factor;
	
	                  factor = bucket_fraction_complement[Peak2D::RT] * bucket_fraction_complement[Peak2D::MZ];
	                  tb( bucket_index[Peak2D::RT], bucket_index[Peak2D::MZ] ) += tq * factor;
	
	                  factor = bucket_fraction_complement[Peak2D::RT] * bucket_fraction[Peak2D::MZ];
	                  tb( bucket_index[Peak2D::RT], bucket_index[Peak2D::MZ] + 1 ) += tq * factor;
	
	                  factor = bucket_fraction[Peak2D::RT] * bucket_fraction_complement[Peak2D::MZ];
	                  tb( bucket_index[Peak2D::RT] + 1, bucket_index[Peak2D::MZ] ) += tq * factor;
	
	                  factor = bucket_fraction[Peak2D::RT] * bucket_fraction[Peak2D::MZ];
	                  tb( bucket_index[Peak2D::RT] + 1, bucket_index[Peak2D::MZ] + 1 ) += tq * factor;
	
	                  ++number_of_considered_element_pairs_for_this_pair_of_buckets;
	
	                  if ( progress_dots &&
	                       ! (number_of_considered_element_pairs_for_this_pair_of_buckets % progress_dots)
	                     )
	                  {
	                    std::cout << 'H' << std::flush;
	                  }
	
	                } // for scene_iter
	              } // for model_iter
	            } // for model_bucket_index_MZ
	          } // for model_bucket_index_RT
	        } // for scene_bucket_index_MZ
	      } // for scene_bucket_index_RT
	
	
	      // Optionally, write debug output as specified in param.
	      if ( getParameters().exists("debug:dump_shift_buckets") )
	      {
	        String dump_filename = getParameters().getValue("debug:dump_shift_buckets");
	        std::ofstream dump_file(dump_filename.c_str());
	        dump_file << "# " << dump_filename << " generated " << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << std::endl;
	        dump_file << "# Shift buckets: xcoord ycoord quality xindex yindex" << std::endl;
	
	        for ( ShiftQualityMatrixType::ConstIterator iter = tb.begin(); iter != tb.end(); ++iter)
	        {
	          std::pair<UInt,UInt> row_col = tb.indexPair(iter-tb.begin());
	          if ( *iter )
	          {
	            dump_file << tbbe_min[Peak2D::RT] + tbs[Peak2D::RT] * row_col.first << ' '
	            << tbbe_min[Peak2D::MZ] + tbs[Peak2D::MZ] * row_col.second << ' '
	            << *iter << ' '
	            << row_col.first << ' '
	            << row_col.second
	            << " #tb" << std::endl ;
	          }
	        }
	        dump_file << "# " << dump_filename << " EOF " << QDateTime::currentDateTime().toString("yyyy-MM-dd hh:mm:ss").toStdString() << std::endl;
	      }
	    }
	
	    /**
	    	@brief Compute the shift.
	
	    	@note shift_buckets_ must have been calculated before.
	    */
	    void computeShift_(TransformationDescription& mapping)
	    {
	      Shift shift;
	
	      // Shorthands ...
	      ShiftQualityMatrixType const & tb = shift_bucket_;
	      DPosition<2> const & tbs = shift_bucket_size_;
	      UInt const (&tbw)[2] = shift_bucket_window_;
	
	      // Find the transformation bucket with highest impact (quality).
	      UInt tb_max_element_index = std::max_element(tb.begin(),tb.end()) - tb.begin();
	      UInt tb_max_indices[2];
	      tb_max_indices[Peak2D::RT] = tb.rowIndex(tb_max_element_index);
	      tb_max_indices[Peak2D::MZ] = tb.colIndex(tb_max_element_index);
	
	      // Compute a weighted average of the shifts nearby the tb_max_element.
	      //Shift result; // initially zero
	
	      DPosition<2> const& tbbe_min = shift_bounding_box_enlarged_.min();
	      int tb_run_indices[2];
	      for ( tb_run_indices[Peak2D::RT]  = std::max ( int (tb_max_indices[Peak2D::RT] - tbw[Peak2D::RT]), 0 );
	            tb_run_indices[Peak2D::RT] <= std::min ( int (tb_max_indices[Peak2D::RT] + tbw[Peak2D::RT]), int (tb.rows()) - 1 );
	            ++tb_run_indices[Peak2D::RT]
	          )
	      {
	        for ( tb_run_indices[Peak2D::MZ]  = std::max ( int (tb_max_indices[Peak2D::MZ] - tbw[Peak2D::MZ]), 0 );
	              tb_run_indices[Peak2D::MZ] <= std::min ( int (tb_max_indices[Peak2D::MZ] + tbw[Peak2D::MZ]), int (tb.cols()) - 1 );
	              ++tb_run_indices[Peak2D::MZ]
	            )
	        {
	          DPosition<2> contribution_position(tbs);
	          for ( UInt dimension = 0; dimension < 2; ++dimension)
	          {
	            contribution_position[dimension] *= tb_run_indices[dimension];
	          }
	          contribution_position += tbbe_min;
	          DoubleReal contribution_quality = tb( tb_run_indices[Peak2D::RT], tb_run_indices[Peak2D::MZ] );
	          shift.quality += contribution_quality;
	          contribution_position *= contribution_quality;
	          shift.position += contribution_position;
	        }
	      }
	      if ( shift.quality != 0 )
	      {
	        shift.position /= -shift.quality ;
	      }
				
				mapping.setName("linear");
	      mapping.setParam("slope",1.0);
	      mapping.setParam("intercept",shift.position[0]);
	    }
	
	    /**
				@brief Compute the shift and similarity for a pair of elements; larger quality values are better.
				
				The returned value should express our confidence that one element might possibly be matched to the other.
				
				Currently this will just calculate the ratio of intensities, either
				"left/right" or "right/left", such that a value between 0 and 1 is returned.
				
				@improvement Take the quality of the elements themselves into account, i.e., how good they fit to their model. (Eva)
	    */
	    Shift shift_( PointType const & left, PointType const & right ) const
	    {
	      Shift shift;
	      shift.position = right.getPosition() - left.getPosition();
	      if ( right.getIntensity() == 0.0 )
	    	{
	    		shift.quality = 0.0;
	      }
	      else
	      {
					DoubleReal result = left.getIntensity() / right.getIntensity();
	      	shift.quality = ( result <= 1.0 ? result : 1. / result );
	      }
	      return shift;
	    }
			
			///Pointer to the model map
			const ConsensusMap* model_map_;
			///Pointer to the scene map
			const ConsensusMap* scene_map_;
			
	    /// Holds the bounding box of all input elements.
	    DBoundingBox<2>  element_map_position_bounding_box_[2];
	
	    /// Holds the enlarged bounding box for all input elements.  It is larger  by about half of a bucket in all directions.
	    DBoundingBox<2>  element_map_position_bounding_box_enlarged_[2];
	
	    /// Holds a bounding box for the input element intensities.
	    DBoundingBox<1> element_map_intensity_bounding_box_[2];
	
	    /// Element indices are stored in theses buckets.
	    ElementBucketMatrixType element_bucket_[2];
	
	    /// Diagonal size of each bucket.
	    DPosition<2> element_bucket_size_;
	
	    /// Shifts are stored (summed up) in these buckets.
	    ShiftQualityMatrixType shift_bucket_;
	
	    /// Holds a bounding box for all possible shift vectors.
	    DBoundingBox<2> shift_bounding_box_;
	
	    /// Holds an enlarged bounding box for all shift vectors.  It is arger by about half of a bucket in all directions.
	    DBoundingBox<2> shift_bounding_box_enlarged_;
	
	    /// Diagonal size of each bucket in shift_bucket_.
	    DPosition<2> shift_bucket_size_;
	
	    /// Number of surrounding buckets of element indices to be considered when computing shifts.
	    UInt element_bucket_window_[2];
	
	    /// Number of surrounding buckets of shift indices to be considered when computing shifts.
	    UInt shift_bucket_window_[2];
  };

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_POSECLUSTERINGSHIFTSUPERIMPOSER_H
