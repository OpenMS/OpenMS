// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiModule.h>
#include <OpenMS/MATH/STATISTICS/AveragePosition.h>

#include <queue>
#include <iostream>
#include <fstream>

namespace OpenMS
{

	/**
	  @brief Simple feature extension algorithm
		
		This algorithm implements the extension phase of the FeatureFinder 
		as described by Groepl et al. (2005)

		We want to determine a region around a seed that is
		provided by the seeder. Initially, this region is
		empty. The boundary of this region is implemented
		using a MutablePriorityQueue which contains only
		the seed at the beginning.

		At each step, we choose a data point from the boundary,
		move it into the region and explore the neigbourhood of
		this point in a cross-wise manner (m/z up, m/z down, rt up
    and rt down). During this exploration we compute the priority
		of all encountered points as a function of the distance from
		the extracted point. If this priority exceeds a threshold,
		we insert the corresponding point into the boundary and proceed.

		We stop the extension phase if all peaks contained in the
    boundary have an intensity lower than a threshold or are too
    distant from the centroid of the feature.

		@image html SimpleExtender.png

		@htmlinclude OpenMS_SimpleExtender.parameters

		@ingroup FeatureFinder
	*/
	template<class PeakType,class FeatureType>
  class SimpleExtender
    : public FeaFiModule<PeakType,FeatureType>,
			public FeatureFinderDefs
  {
  	public:
		typedef FeaFiModule<PeakType,FeatureType> Base;

		/// Intensity of a data point
  	typedef typename Base::IntensityType IntensityType;
		/// Coordinates of a point (m/z and rt)
  	typedef typename Base::CoordinateType CoordinateType;
		/// Priority of a point (see below)
  	typedef DoubleReal ProbabilityType;

  	/// Constructor
    SimpleExtender(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff)
		: Base(map,features,ff),
			last_pos_extracted_()
		{
			this->setName("SimpleExtender");
      
      this->defaults_.setValue("dist_mz_up",6.0,"Maximum high m/z distance of peak in the region/boundary from the seed.");
      this->defaults_.setMinFloat("dist_mz_up",0.0);
			this->defaults_.setValue("dist_mz_down",2.0,"Maximum low m/z distance of peak in the region/boundary from the seed.");
      this->defaults_.setMinFloat("dist_mz_down",0.0);
			this->defaults_.setValue("dist_rt_up",5.0,"Maximum high RT distance of peak in the region/boundary from the seed.");
      this->defaults_.setMinFloat("dist_rt_up",0.0);
			this->defaults_.setValue("dist_rt_down",5.0,"Maximum low RT distance of peak in the region/boundary from the seed.");
      this->defaults_.setMinFloat("dist_rt_down",0.0);

			// priority check is per default switched off
			// these values were used for the Myoglobin quantification project
			// DON'T REMOVE THIS
			this->defaults_.setValue("priority_thr",-0.1,"Minimum priority for data points to be included into the boundary of the feature (default 0.0). The priority of a data point is a function of its intensity and its distance to the last point included into the feature region. Setting this threshold to zero or a very small value is usually a good idea.", StringList::create("advanced"));
     
			this->defaults_.setValue("intensity_factor",0.03,"Influences for intensity (ion count) threshold in the feature extension. We include only raw data points into this region if their intensity is larger than [intensity_factor * (intensity of the seed)].");
      this->defaults_.setMinFloat("intensity_factor",0.0);
      this->defaults_.setMaxFloat("intensity_factor",1.0);
      
			this->defaultsToParam_();
		}

    /// destructor
    virtual ~SimpleExtender()
		{
		}

    /// return next seed
    void extend(const ChargedIndexSet& seed_region, ChargedIndexSet& result_region)
		{
			// empty region and boundary datastructures
			result_region.clear();
			priorities_.clear();
			running_avg_.clear();
			boundary_ = std::priority_queue< IndexWithPriority, std::vector<IndexWithPriority>, typename IndexWithPriority::PriorityLess>();

#ifdef DEBUG_FEATUREFINDER
			std::vector<IndexPair> debug_vector;
#endif
			// find maximum of region (seed)
			CoordinateType max_intensity = 0.0;
			IndexPair seed;

			for (IndexSet::const_iterator citer = seed_region.begin(); citer != seed_region.end(); ++citer)
			{
				if (this->getPeakIntensity(*citer) > max_intensity)
				{
					seed = *citer;
					max_intensity = this->getPeakIntensity(seed);
				}
			}

			// remember last extracted point (in this case the seed !)
			last_pos_extracted_[Peak2D::RT] = this->getPeakRt(seed);
			last_pos_extracted_[Peak2D::MZ] = this->getPeakMz(seed);

			// Add peaks received from seeder directly to boundary
			for (IndexSet::const_iterator citer = seed_region.begin(); citer != seed_region.end(); ++citer)
			{
				ProbabilityType priority = computePeakPriority_(*citer);
				priorities_[*citer] = priority;
				boundary_.push(IndexWithPriority(*citer,priority));
			}
			// pass on charge information
			result_region.charge = seed_region.charge;

			// re-compute intensity threshold
			intensity_threshold_ = (DoubleReal)(this->param_).getValue("intensity_factor") * this->getPeakIntensity(seed);

#ifdef DEBUG_FEATUREFINDER
			std::cout << "\n";
			std::cout << "Extending from " << this->getPeakRt(seed) << "/" << this->getPeakMz(seed) << std::endl;
			std::cout << "Intensity of seed " << this->getPeakIntensity(seed);
			std::cout << " (" << seed.first << "/" << seed.second << ")" << std::endl;
			std::cout << "Intensity_threshold: " << intensity_threshold_ << std::endl;
#endif

			while (!boundary_.empty())
			{
				// remove peak with highest priority
				const IndexPair  current_index = boundary_.top().index;
				boundary_.pop();

				// 	check for corrupt index
				OPENMS_PRECONDITION(current_index.first<(*this->map_).size(), "Scan index outside of map!");
				OPENMS_PRECONDITION(current_index.second<(*this->map_)[current_index.first].size(), "Peak index outside of scan!");

				// remember last extracted peak
				last_pos_extracted_[Peak2D::RT] = this->getPeakRt(current_index);
				last_pos_extracted_[Peak2D::MZ] = this->getPeakMz(current_index);

				// Now we explore the neighbourhood of the current peak. Points in this area are included
				// into the boundary if their intensity is not too low and they are not too
				// far away from the seed.
				// Add position to the current average of positions weighted by intensity
				running_avg_.add(last_pos_extracted_,this->getPeakIntensity(current_index));

				// explore neighbourhood of current peak
				moveMzUp_(current_index);
				moveMzDown_(current_index);
				moveRtUp_(current_index);
				moveRtDown_(current_index);

				// set peak flags and add to boundary
				this->ff_->getPeakFlag(current_index) = USED;
#ifdef DEBUG_FEATUREFINDER
				debug_vector.push_back(current_index);
#endif
				result_region.insert(current_index);

			} // end of while ( !boundary_.empty() )

#ifdef DEBUG_FEATUREFINDER
			std::cout << "Feature region size: " << result_region.size() << std::endl;
#endif

#ifdef DEBUG_FEATUREFINDER
			static UInt number=1;
			writeDebugFile_(debug_vector,number++);
			debug_vector.clear();
#endif

			return;
		} // end of extend

    /**
     @brief A helper structure to sort indizes by their priority.

     This structure is used to keep track of the boundary of a
     feature. After a peak is found during the extension phase,
     we compute its priority (which is dependant on its distance from
     the point that was the last to be extracted from the boundary
     and its intensity). If this priority is large enough, we include
     the point into the boundary. The boundary (which is implemented
     as mutable priority queue) sorts the peaks by this priority.

    */
  	struct IndexWithPriority
  	{
  		IndexWithPriority(const FeatureFinderDefs::IndexPair& i, DoubleReal p) 
				: index(i), 
					priority(p)
  		{
  		}

  		IndexPair index;
  		ProbabilityType priority;

			///Compares two indizes by priority.
  		struct PriorityLess
  		{
  			inline bool operator() (const IndexWithPriority& x, const IndexWithPriority& y) const
				{
    			return x.priority < y.priority;
				}
			};
  	};

  protected:
  
  	virtual void updateMembers_()
		{
			dist_mz_up_ = this->param_.getValue("dist_mz_up");
			dist_mz_down_ = this->param_.getValue("dist_mz_down");
			dist_rt_up_ = this->param_.getValue("dist_rt_up");
			dist_rt_down_ = this->param_.getValue("dist_rt_down");
			priority_threshold_ = this->param_.getValue("priority_thr");
		}

		/// write DTA2D debug file for the feature with index @p nr_feat
  	void writeDebugFile_(const std::vector<IndexPair>& peaks, UInt nr_feat)
		{
			String filename = String(nr_feat).fillLeft('0',4) + "_Extension.dta2d";
			std::ofstream file(filename.c_str());
			for (Size i=0; i<peaks.size(); ++i)
			{
				file << this->getPeakRt(peaks[i]) << " " << this->getPeakMz(peaks[i]) << " " << peaks.size()-i << std::endl;
			}
			file.close();
		}

  	/// Checks if the current peak is too far from the centroid
  	bool isTooFarFromCentroid_(const IndexPair& index)
		{
			//Corrupt index
			OPENMS_PRECONDITION(index.first < (*this->map_).size(), "Scan index outside of map!");
			OPENMS_PRECONDITION(index.second < (*this->map_)[index.first].size() , "Peak index outside of scan!");

			const DPosition<2>& curr_mean = running_avg_.getPosition();

			if ( this->getPeakMz(index) > curr_mean[Peak2D::MZ] + dist_mz_up_   ||
					 this->getPeakMz(index) < curr_mean[Peak2D::MZ] - dist_mz_down_ ||
					 this->getPeakRt(index) > curr_mean[Peak2D::RT] + dist_rt_up_   ||
					 this->getPeakRt(index) < curr_mean[Peak2D::RT] - dist_rt_down_ )
			{
				//too far
				return true;
			}

			//close enough
			return false;
		}

   	/// Extends the seed into positive m/z direction
  	void moveMzUp_(const IndexPair& index)
		{
			try
			{
				IndexPair tmp = index;
				while (true)
				{
					this->getNextMz(tmp);
					if (isTooFarFromCentroid_(tmp)) break;
					checkNeighbour_(tmp);
				}
			}
			catch(NoSuccessor)
			{
			}
		}

  	/// Extends the seed into negative m/z direction
  	void moveMzDown_(const IndexPair& index)
		{
			try
			{
				IndexPair tmp = index;
				while (true)
				{
					this->getPrevMz(tmp);
					if (isTooFarFromCentroid_(tmp))	break;
					checkNeighbour_(tmp);
				}
			}
			catch(NoSuccessor)
			{
			}
		}

  	/// Extension into positive rt dimension
  	void moveRtUp_(const IndexPair& index)
		{
			try
			{
				IndexPair tmp = index;

				while (true)
				{
					this->getNextRt(tmp);
					if (isTooFarFromCentroid_(tmp)) break;
					checkNeighbour_(tmp);
				}
			}
			catch(NoSuccessor)
			{
			}
		}

  	/// Extends the seed into negative retention time direction
  	void moveRtDown_(const IndexPair& index)
		{
			try
			{
				IndexPair tmp = index;
				while (true)
				{
					this->getPrevRt(tmp);
					if (isTooFarFromCentroid_(tmp)) break;
					checkNeighbour_(tmp);
				}
			}
			catch(NoSuccessor)
			{
			}
		}

  	/// Computes the priority of a peak as function of intensity and distance from seed.
  	ProbabilityType computePeakPriority_(const IndexPair& index)
		{
     return (*this->map_)[index.first][index.second].getIntensity();
    }

  	/// Checks the neighbours of the current for insertion into the boundary.
  	void checkNeighbour_(const IndexPair& index)
		{
			//Corrupt index
			OPENMS_PRECONDITION(index.first<(*this->map_).size(), "Scan index outside of map!");
			OPENMS_PRECONDITION(index.second<(*this->map_)[index.first].size(), "Peak index outside of scan!");

			// skip this point if its intensity is too low
			if (this->getPeakIntensity(index) <= intensity_threshold_)
			{
			 return;
			}
			if ( this->ff_->getPeakFlag(index) == UNUSED)
			{
				DoubleReal pr_new = computePeakPriority_(index);

				if (pr_new > priority_threshold_)
				{
					//std::map<IndexPair, DoubleReal>::iterator piter = priorities_.find(index);
					this->ff_->getPeakFlag(index) = USED;
					priorities_[index] = pr_new;
					boundary_.push(IndexWithPriority(index,pr_new));
				}
  		}
		}

  	/// keeps an running average of the peak coordinates weighted by the intensities
  	Math::AveragePosition<2> running_avg_;

  	/// Keeps track of peaks already included in the boundary (value is priority of peak)
  	std::map<IndexPair, ProbabilityType> priorities_;

  	/// Position of last peak extracted from the boundary (used to compute the priority of neighbouring peaks)
  	DPosition<2> last_pos_extracted_;

  	/// Represents the boundary of a feature
  	std::priority_queue< IndexWithPriority, std::vector < IndexWithPriority > , typename IndexWithPriority::PriorityLess > boundary_;

		/// Mininum intensity of a boundary point. Calculated from 'intensity_factor' and the seed intensity
		IntensityType intensity_threshold_;

		/// Maximum distance to seed in positive m/z
		CoordinateType dist_mz_up_;
		/// Maximum distance to seed in negative m/z
		CoordinateType dist_mz_down_;
		/// Maximum distance to seed in positive retention time
		CoordinateType dist_rt_up_;
		/// Maximum distance to seed in negative retention time
		CoordinateType dist_rt_down_;

		/// Minium priority for points in the feature region (priority is function of intensity and distance to seed)
		ProbabilityType priority_threshold_;
		
    /// charged index set
		ChargedIndexSet region_;

		private:
			/// Not implemented
			SimpleExtender();
			/// Not implemented
			SimpleExtender& operator=(const SimpleExtender&);
			/// Not implemented
			SimpleExtender(const SimpleExtender&);

  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_SIMPLEEXTENDER_H
