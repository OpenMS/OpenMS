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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SimpleExtender.h>

namespace OpenMS
{

	SimpleExtender::SimpleExtender() : BaseExtender(),
																		 first_seed_seen_(false), intensity_threshold_(0),
																		 last_extracted_(0), nr_peaks_seen_(0)
	{
    name_ = SimpleExtender::getName();
    defaults_.setValue("tolerance_rt",2.0f);
    defaults_.setValue("tolerance_mz",0.5f);
    defaults_.setValue("dist_mz_up",6.0f);
    defaults_.setValue("dist_mz_down",2.0f);
    defaults_.setValue("dist_rt_up",5.0f);
    defaults_.setValue("dist_rt_down",5.0f);
    defaults_.setValue("priority_thr",0.01f);
    defaults_.setValue("intensity_factor",0.03f);

    param_ = defaults_;
	}

	SimpleExtender::~SimpleExtender()
	{}

	const IndexSet& SimpleExtender::extend(const UnsignedInt seed_index)
	{
    if (!first_seed_seen_)
    {

			float tol_rt = param_.getValue("tolerance_rt");
			float tol_mz = param_.getValue("tolerance_mz");
			intensity_factor_ = param_.getValue("intensity_factor");

			dist_mz_up_   = param_.getValue("dist_mz_up");
			dist_mz_down_ = param_.getValue("dist_mz_down");
			dist_rt_up_   = param_.getValue("dist_rt_up");
			dist_rt_down_ = param_.getValue("dist_rt_down");

			priority_threshold_ = param_.getValue("priority_thr");

			// initialise distributions to interpolate the priority
			score_distribution_rt_.getData().push_back( 1. );
			score_distribution_rt_.setScale(tol_rt);
			score_distribution_rt_.setOffset(0);

			score_distribution_mz_.getData().push_back(1.);
			score_distribution_mz_.setScale(tol_mz);
			score_distribution_mz_.setOffset(0);

			// reserve enough space in mutable queue
			//unsigned int peak_nr = traits_->getNumberOfPeaks();
			//boundary_.reserve(peak_nr);

			first_seed_seen_ = true;

    }

    // empty region and boundary datastructures
    region_.clear();
    //boundary_.clear(); // works only with mutable queue
    while (boundary_.size() > 0)
			boundary_.pop();
    priorities_.clear();
    running_avg_.clear();

    nr_peaks_seen_       = 0;
    intensity_threshold_ = 0;

    // remember last peak to be extracted from the boundary
    // in this case, the seed
    last_extracted_ = seed_index;

    float prior = computePeakPriority_(seed_index);
    IndexWithPriority seed(seed_index,prior);

    // at the beginning, the boundary contains only the seed
    boundary_.push(seed);
    priorities_[seed_index] = prior;

    while (!boundary_.empty())
    {
			// remove peak with highest priority
			IndexWithPriority const index_priority = boundary_.top();
			boundary_.pop();

			nr_peaks_seen_++;

			UnsignedInt const current_index       = index_priority.index;
			IntensityType const current_intensity = traits_->getPeakIntensity(current_index);

			// remember last peak to be extracted from the boundary
			last_extracted_ = current_index;

			// we set the intensity threshold for raw data points to be
			// included into the feature region to be a fraction of the
			// intensity of the fifth largest peak
			if (nr_peaks_seen_ == 5)
				intensity_threshold_ = intensity_factor_ * current_intensity;

			if (current_intensity <  intensity_threshold_)
				continue;

			// Now we explore the neighbourhood of the current peak. Peaks in this area are included
			// into the boundary if their intensity is not too low and they are not too
			// far away from the seed.

			// get position of current peak
			FeaFiTraits::PositionType const curr_pos = traits_->getPeak(current_index).getPosition();
			// and add it to the current average of positions weighted by the intensity
			running_avg_.add(curr_pos,current_intensity);

			// explore neighbourhood of current peak
			moveMzUp_(current_index);
			moveMzDown_(current_index);
			moveRtUp_(current_index);
			moveRtDown_(current_index);

			// set flag
			if (traits_->getPeakFlag(current_index) != FeaFiTraits::SEED)
				traits_->getPeakFlag(current_index) = FeaFiTraits::INSIDE_FEATURE;

			region_.add(current_index);

    } // end of while ( !boundary_.empty() )

    std::cout << "Feature region size: " << region_.size() << std::endl;

    region_.sort();
    return region_;

	} // end of extend(Unsigned int seed_index)

	/**
	 * \brief Checks whether the current peak is to far away from the seed
	 * 
	 * */
	bool SimpleExtender::isTooFarFromCentroid_(UnsignedInt current_peak)
	{

    CoordinateType const current_mz   = traits_->getPeakMz(current_peak);
    CoordinateType const current_rt   = traits_->getPeakRt(current_peak);

    FeaFiTraits::PositionType const curr_mean = running_avg_.getPosition();

    float dist_mz_up = param_.getValue("dist_mz_up");
    float dist_mz_down = param_.getValue("dist_mz_down");
    float dist_rt_up = param_.getValue("dist_rt_up");
    float dist_rt_down = param_.getValue("dist_rt_down");

    if ( current_mz > curr_mean[MZ] + dist_mz_up   ||
				 current_mz < curr_mean[MZ] - dist_mz_down ||
				 current_rt > curr_mean[RT] + dist_rt_up   ||
				 current_rt < curr_mean[RT] - dist_rt_down )
    {
			return true;
    }
    return false;
	}

	void SimpleExtender::moveMzUp_(UnsignedInt current_index)
	{

    UnsignedInt current_scan = traits_->getPeakScanNr(current_index);

    try // moving up in m/z direction
    {
			while (true)
			{
				current_index = traits_->getNextMz(current_index); // take next peak

				// ???? TODO (Clemens)
				// Various improvements are possible here:
				// The following triggers a binary search in scan_index_.
				// This could be avoided by looking at RT itself instead of ScanNr.
				// Also, isTooFarFromCentroid_ need not consider RT here, since it has not changed.

				// stop if we've left the current scan
				if (current_scan != traits_->getPeakScanNr(current_index)
						|| isTooFarFromCentroid_(current_index) )
					break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(current_index);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {}
	}

	void SimpleExtender::moveMzDown_(UnsignedInt current_index)
	{

    UnsignedInt current_scan = traits_->getPeakScanNr(current_index);

    try // moving down in m/z direction
    {
			while (true)
			{
				current_index    = traits_->getPrevMz(current_index);	// take next peak

				// ???? TODO (Clemens) 
				// see remarks in moveMzUp_()

				// stop if we've left the current scan
				if (current_scan != traits_->getPeakScanNr(current_index)
						|| isTooFarFromCentroid_(current_index))
					break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(current_index);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {}
	}

	void SimpleExtender::moveRtUp_(UnsignedInt current_index)
	{

    try // moving up in retention time
    {
			while (true)
			{

				// ???? TODO (Clemens)
				// Again, binary search for rt is not necessary.  We could maintain that information.
				// Also, binary search within scan should look for mz of staring point, not the last mz.
				// (Note: isTooFarFromCentroid_() _does_ need to consider both rt and mz here, since mz might have changed.)

				current_index = traits_->getNextRt(current_index); // take next peak

				if (isTooFarFromCentroid_(current_index))
					break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(current_index);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {}
	}

	void SimpleExtender::moveRtDown_(UnsignedInt current_index)
	{

		// ???? TODO (Clemens) 
		// see remarks in moveRtUp_()


    try // moving up in retention time
    {

			while (true)
			{
				current_index = traits_->getPrevRt(current_index); // take next peak

				if (isTooFarFromCentroid_(current_index))
					break;

				// check this neighbour for insertion into the boundary
				checkNeighbour_(current_index);

			} // end of while (true)

    }
    catch(NoSuccessor)
    {}
	}

	double SimpleExtender::computePeakPriority_(UnsignedInt current_peak)
	{

    IntensityType curr_intens = traits_->getPeakIntensity(current_peak);

    CoordinateType curr_mz = traits_->getPeakMz(current_peak);
    CoordinateType last_mz = traits_->getPeakMz(last_extracted_);

    CoordinateType curr_rt = traits_->getPeakRt(current_peak);
    CoordinateType last_rt = traits_->getPeakRt(last_extracted_);

    return curr_intens *
			score_distribution_rt_.value(curr_rt-last_rt) *
			score_distribution_mz_.value(curr_mz-last_mz);

	}

	void SimpleExtender::checkNeighbour_(UnsignedInt current_index)
	{
    // we don't care about points with intensity zero
    if (traits_->getPeakIntensity(current_index) == 0)
			return;

    if (traits_->getPeakFlag(current_index)== FeaFiTraits::UNUSED
				/*&& !isTooFarFromCentroid_(current_index)*/ )
    {
			std::map<UnsignedInt, double>::iterator piter = priorities_.find(current_index);
			double pr_new = computePeakPriority_(current_index);

			if (piter == priorities_.end()) // not yet in boundary
			{
				if (pr_new > priority_threshold_) // check if priority larger than threshold
				{
					priorities_[current_index] = pr_new;
					IndexWithPriority peak(current_index,pr_new);
					boundary_.push(peak);	// add to boundary
				}
			}
			/*else // already in boundary
        {
				double pr_curr = piter->second;
				if (pr_new > pr_curr) // update priority only if larger than old one.
				{
				priorities_[current_index] = pr_new;
				IndexWithPriority peak(current_index,pr_new);
				boundary_.update(peak);	 // update boundary
				} 
        } // end if (piter == priorities_.end())
			*/
    } // end if traits_->...
	}


} // end of class SimpleExtender
