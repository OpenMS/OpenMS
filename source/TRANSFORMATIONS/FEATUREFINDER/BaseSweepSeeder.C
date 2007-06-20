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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSweepSeeder.h>

#include<OpenMS/SYSTEM/StopWatch.h>

using namespace std;

namespace OpenMS
{

BaseSweepSeeder::BaseSweepSeeder()
        : BaseSeeder(),
				is_initialized_(false),
				mass_tolerance_alignment_(0),
				scans_to_sumup_(0),
				mass_tolerance_cluster_(0)
{
		// number of scans used for alignment
		defaults_.setValue("scans_to_sumup",5);
		// mass tolerance during scan alignment
		defaults_.setValue("mass_tolerance_alignment", 0.1);
		
			// minimum number of scans per isotopic cluster
		defaults_.setValue("min_number_scans",5);
			// maximum number of scans per isotopic cluster
		defaults_.setValue("max_number_scans",300);
		// minimum number of peaks per cluster
		defaults_.setValue("min_number_peaks",20);		
		
		// mass tolerance for point cluster construction
		defaults_.setValue("mass_tolerance_cluster",1.2);
		// rt tolerance for cluster construction (given in number of scans)
		defaults_.setValue("rt_tolerance_cluster",2);		
		
		// max distance in rt for merged peak cluster (given in # scans)
		defaults_.setValue("max_rt_dist_merging",40.0);
		// max distance in mz for merged peak cluster 
		defaults_.setValue("max_mz_dist_merging",1.5);
}

BaseSweepSeeder::BaseSweepSeeder(const BaseSweepSeeder& source) : BaseSeeder(source) {}

BaseSweepSeeder::~BaseSweepSeeder() {}

BaseSweepSeeder& BaseSweepSeeder::operator = (const BaseSweepSeeder& source)
{
    if (&source == this)
        return *this;

    BaseSeeder::operator = (source);

    return *this;
}

FeaFiModule::ChargedIndexSet BaseSweepSeeder::nextSeed() throw (NoSuccessor)
{
		if (!is_initialized_)
		{
			StopWatch w2;
			w2.start();
			sweep_();		// sweep across map and scan for pattern 
			curr_region_  = iso_map_.begin();
			is_initialized_ = true;
			w2.stop();
			cout << "Seeding took : " << w2.getClockTime() << " [s]" << endl;
		}
		
		if ( curr_region_ == iso_map_.end() || iso_map_.size() == 0 )
		{
			throw NoSuccessor(__FILE__, __LINE__,__PRETTY_FUNCTION__, make_pair(0,0));
		}
		
		cout << "Retrieving next region with charge " << (*curr_region_).second.peaks_.charge_; 
		cout << " and size " << (*curr_region_).second.peaks_.size() << endl;
			
		return (curr_region_++)->second.peaks_;
}

void BaseSweepSeeder::updateMembers_()
{
	// params for scan alignment
	mass_tolerance_alignment_ = param_.getValue("mass_tolerance_alignment");
	scans_to_sumup_               = param_.getValue("scans_to_sumup");

	// sweepline params
	mass_tolerance_cluster_      = param_.getValue("mass_tolerance_cluster");
	rt_tolerance_cluster_            = (UInt) param_.getValue("rt_tolerance_cluster");
	
	// cluster merging params
	max_rt_dist_merging_          = param_.getValue("max_rt_dist_merging"); 
	max_mz_dist_merging_       = param_.getValue("max_mz_dist_merging");
}

void BaseSweepSeeder::sweep_()
{
		// progress logger
		traits_->startProgress(0, traits_->getData().size() , "FeatureFinder");
		
		// copy current scan. 
		// This is necessary as the peak intensities have to be modified in the sumUp_ method
 		//SpectrumType current_scan  = traits_->getData()[0];
				
		for (UInt currscan_index = 0; currscan_index < traits_->getData().size(); ++currscan_index)
		{		
			traits_->setProgress(currscan_index);
			
			SpectrumType current_scan  = traits_->getData()[ currscan_index ];		
			
			cout << "---------------------------------------------------------------------------" << endl;
			cout << "Processing scan " << (currscan_index + 1) << " of " << traits_->getData().size() << endl;
			cout << "Retention time: " << current_scan.getRT() << endl;
			
			#ifdef DEBUG_FEATUREFINDER
			// write debug output
			String fname = String("scan_") + current_scan.getRT();;
			ofstream out( fname.c_str() );
			for(UInt k = 0; k<current_scan.size();++k)
			{
				out << current_scan[k].getMZ() << " " << current_scan[k].getIntensity() << endl;
			}
			out.close();
			#endif
			
			StopWatch w;
			w.start();			
			sumUp_(current_scan,currscan_index);
			w.stop();
			cout << "Alignment took " << w.getClockTime() << " [s]. " << endl;
			w.reset();
// 			if (currscan_index == 0)
// 			{
// 				// align and sum up first scan(s)
// 				sumUp_(current_scan,currscan_index);
// 			}
// 			else
// 			{
// 				// add next one
// 				addNextScan_(current_scan,currscan_index);			
// 				// substract last scan
// 				substractLastScan_(current_scan,currscan_index);
// 			}
						
			#ifdef DEBUG_FEATUREFINDER
			// write debug output
			fname = String("scan_aligned_") + current_scan.getRT();;
			out.open( fname.c_str() );
			for(UInt k = 0; k<current_scan.size();++k)
			{
				out << current_scan[k].getMZ() << " " << current_scan[k].getIntensity() << endl;
			}
			out.close();
			#endif
						
			// detect isotopic pattern...
			w.start();		
			ScoredMZVector iso_curr_scan = detectIsotopicPattern_(current_scan );
			w.stop();
			cout << "Isotopic pattern detection took " << w.getClockTime() << " [s]. " << endl;
			w.reset();	
			
			for (ScoredMZVector::const_iterator citer = iso_curr_scan.begin();
						citer != iso_curr_scan.end();
						++citer)
			{
 				traits_->getPeakFlag( make_pair( currscan_index, citer->first ) ) = FeaFiTraits::USED;			
			}
			
			// for each m/z position with score: 
			// => check for cluster at similar m/z in previous scans
			// => if (matching cluster found) extend
			for (ScoredMZVector::const_iterator citer = iso_curr_scan.begin();
						citer != iso_curr_scan.end();
						++citer)
			{
				// check if we have another cluster close by
				TableIteratorType entry_to_insert = checkInPreviousScans_(*citer,currscan_index);
				
				#ifdef DEBUG_FEATUREFINDER
			  cout << "Cluster at (" << traits_->getPeakRt( make_pair( currscan_index, citer->first ) );
				cout << " / " << traits_->getPeakMz( make_pair( currscan_index, citer->first ) ) << ")" << endl;
				cout << "Charge estimate: " << citer->second.first << " score " << citer->second.second << endl;
				cout << "Extending...." << endl;
				#endif
							
				// store charge estimate and score in this scan
				entry_to_insert->second.scored_charges_.push_back( citer->second );
				// store seed
				entry_to_insert->second.peaks_.insert( make_pair(currscan_index,citer->first) );
				
				// Problem: how to perform the extension (i.e. handle scans with weak signals)
				// Solution: use very low threshold in cwt, return many possible positions, use high threshold if exact match in previous
				// scan, lower threshold if previous scan is further away.
				UInt this_peak                =  citer->first;
				CoordinateType start_mz = traits_->getPeakMz( make_pair(currscan_index,this_peak) );
				CoordinateType mz_dist  = 0;
						
				// walk to the left (for at most 3 Th)
				while (mz_dist < 3.0 && this_peak >= 1)
				{
					mz_dist = ( start_mz - traits_->getPeakMz( make_pair(currscan_index,this_peak) ) );
					
					if ( traits_->getPeakFlag( make_pair(currscan_index,this_peak) ) == FeaFiTraits::UNUSED )
					{
						entry_to_insert->second.peaks_.insert( make_pair(currscan_index,this_peak) );
						traits_->getPeakFlag( make_pair(currscan_index,this_peak) ) = FeaFiTraits::USED;
					}
					--this_peak;
				}
			
				// reset
				this_peak =  (citer->first+1);
				mz_dist   = ( traits_->getPeakMz( make_pair(currscan_index,this_peak) )  - start_mz );
					
				// and to the right (we walk for at most 3 Th)
				while (mz_dist < 2.0 && this_peak < current_scan.size() )
				{
					if ( traits_->getPeakFlag( make_pair(currscan_index,this_peak) )  == FeaFiTraits::UNUSED )
					{
						entry_to_insert->second.peaks_.insert( make_pair(currscan_index,this_peak) );
						traits_->getPeakFlag( make_pair(currscan_index,this_peak) ) = FeaFiTraits::USED;
					}					
					mz_dist = ( traits_->getPeakMz( make_pair(currscan_index,++this_peak) )  - start_mz );
				}
							
			}
					
		}	// end loop for all scans
		
		#ifdef DEBUG_FEATUREFINDER 
		cout << "-----------------------------------------------------------" << endl;
		cout << "List of seeding regions: " << endl;
		for (TableConstIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
		{
			cout << "m/z " << iter->first << " charge: " << iter->second.peaks_.charge_ ;
			cout << " first scan " << iter->second.first_scan_ << " last scan " << iter->second.last_scan_ << endl;
			
			for (vector<UInt>::const_iterator citer = 	iter->second.scans_.begin(); 
						citer != iter->second.scans_.end();
						++citer)
			{
				cout << "# scan : " << *citer << " (";
				IDX tmp;
				tmp.first = *citer;
				cout << traits_->getPeakRt(tmp) << ")" << endl;				
			}
		
		}
		#endif
		

		// fliter hash entries (by number of scans and number of points in the cluster)
		filterHash_();		
		
		// determine most likely charge state by majority voting
		voteForCharge_();
		
		// debug output of all seeding regions with charge
		#ifdef DEBUG_FEATUREFINDER 
		cout << "-----------------------------------------------------------" << endl;
		cout << "List of seeding regions: " << endl;
		for (TableConstIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
		{
			cout << "m/z " << iter->first << " charge: " << iter->second.peaks_.charge_ ;
			cout << " first scan " << iter->second.first_scan_ << " last scan " << iter->second.last_scan_ << endl;
			
			for (vector<UInt>::const_iterator citer = 	iter->second.scans_.begin(); 
						citer != iter->second.scans_.end();
						++citer)
			{
				cout << "# scan : " << *citer << " (";
				IDX tmp;
				tmp.first = *citer;
				cout << traits_->getPeakRt(tmp) << ")" << endl;				
			}
		
		}
		#endif
		
}

void BaseSweepSeeder::computeBorders_(TableIteratorType& entry)
{
	if (entry->second.scans_.size() == 0) return;		
	
	sort(entry->second.scans_.begin(),entry->second.scans_.end());
	
	entry->second.first_scan_ = *entry->second.scans_.begin();
	entry->second.last_scan_ = *(entry->second.scans_.end() - 1);
}

void BaseSweepSeeder::deleteHashEntries_(std::vector<TableIteratorType>& entries)
{
	//cout << "Deleting " << entries.size() << " entries. " << endl;
	for (UInt i=0; i<entries.size();++i)
	{
		iso_map_.erase(entries[i]);
	}

}

void BaseSweepSeeder::filterForOverlaps_()
{
		if (iso_map_.size() == 0) return; // nothing to do

		vector<bool> seen(iso_map_.size(),false);		
		
		vector<TableIteratorType> entries_to_delete;
		vector<Int> indizes;

		UInt counter = 0;
		
		for (TableIteratorType iter = iso_map_.begin(); counter < (iso_map_.size() - 1); ++iter, ++counter)
		{			
			TableIteratorType tmp_iter = iter;			
			++tmp_iter;
		
			CoordinateType mz_dist = tmp_iter->first - iter->first;
			UInt rt_dist                    = 0;

			if (seen.at(counter)) 
			{
				continue; // we've seen that one already
			}
			UInt tmpc = counter + 1;
			
			while (  tmp_iter !=iso_map_.end() && mz_dist < max_mz_dist_merging_ )
			{			
				bool rt_overlap = false;
				
				// test if cluster overlap
				if (	(iter->second.first_scan_ >= tmp_iter->second.first_scan_ && // first case
				       iter->second.first_scan_ <= tmp_iter->second.last_scan_ ) ||
							(iter->second.last_scan_ >= tmp_iter->second.first_scan_ && // second case
				       iter->second.last_scan_ <= tmp_iter->second.first_scan_ ) ||
							(tmp_iter->second.first_scan_ >= iter->second.first_scan_ && // third case
				       tmp_iter->second.first_scan_ <= iter->second.last_scan_ ) ||
							(tmp_iter->second.last_scan_ >= iter->second.first_scan_ && // fourth case
				       tmp_iter->second.last_scan_ <= iter->second.first_scan_ ) )
				{
					rt_overlap = true;				
				}
				       
				rt_dist   = tmp_iter->second.first_scan_ - iter->second.last_scan_;
				mz_dist = tmp_iter->first - iter->first;
				
				//if (rt_dist < max_rt_dist_merging_) rt_overlap = true;
						
				// we merge only features with the same charge, overlap in rt and if we haven't seen them yet.
				if ( tmp_iter->second.peaks_.charge_ == iter->second.peaks_.charge_ 
				     && rt_overlap
						 && !seen.at(tmpc)) 
				{
					// merging features...
					// copy scans
					for (std::vector<UInt>::const_iterator scan_iter = tmp_iter->second.scans_.begin(); 
								scan_iter != tmp_iter->second.scans_.end();
								++scan_iter) 
					{
						if ( *(iter->second.scans_.end() - 1 ) != *scan_iter)	// scan already contained ??
						{
							iter->second.scans_.push_back(*scan_iter);
						}
					}
					
					// copy peaks
					for (std::set<IDX>::const_iterator set_iter = tmp_iter->second.peaks_.begin();
								set_iter != tmp_iter->second.peaks_.end();
								++set_iter)
					{
						iter->second.peaks_.insert(*set_iter);	
					}
					
					// copy charge estimates for each scan
					for (std::vector< ScoredChargeType >::const_iterator sc_iter =  tmp_iter->second.scored_charges_.begin();
								sc_iter !=  tmp_iter->second.scored_charges_.end();
								++sc_iter)
					{
						iter->second.scored_charges_.push_back(*sc_iter);					
					}
					
					// recompute median
					computeBorders_(iter);		
					
					// mark this cluster as deleted
					entries_to_delete.push_back(tmp_iter);
					indizes.push_back(tmpc);
					// and seen
					seen.at(tmpc) = true;
				} // end of if (...)
				
				++tmp_iter;
				++tmpc;
			} // end of while
						
		} // end of for all table entries
		
		deleteHashEntries_(entries_to_delete);
			
} // end of filterHashForOverlaps_(...)


void BaseSweepSeeder::filterForSize_()
{
		// filter for number of scans / significance
		vector<TableIteratorType> entries_to_delete;
	
		UInt min_number_scans = param_.getValue("min_number_scans");
		UInt max_number_scans = param_.getValue("max_number_scans");
		
		UInt min_number_peaks = param_.getValue("min_number_peaks");
				
		// Filter point cluster
		for (TableIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
		{				
			std::vector<UInt>::iterator new_end = std::unique(iter->second.scans_.begin(),iter->second.scans_.end());
			iter->second.scans_.erase(new_end,iter->second.scans_.end());
		
			if (iter->second.scans_.size() < min_number_scans || 
					iter->second.scans_.size() > max_number_scans || 
			    iter->second.peaks_.size() < min_number_peaks)
			{
				entries_to_delete.push_back(iter);
			}
		}	
		
		deleteHashEntries_(entries_to_delete);		
}

void BaseSweepSeeder::filterForSignificance_()
{
	vector<TableIteratorType> entries_to_delete;
		
	ProbabilityType alpha = 0.2;
	//cout << "Filtering for significance with alpha: " << (alpha/iso_map_.size() ) << endl;	

	for (TableIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
	{				
		if (iter->second.peaks_.max_charge_score_ > (alpha/iso_map_.size() ) )
		{
			entries_to_delete.push_back(iter);	
		}
	}
	
	deleteHashEntries_(entries_to_delete);	
}


void BaseSweepSeeder::filterHash_()
{
		cout << iso_map_.size() << " isotopic clusters were found." << endl;
					
		for (TableIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
		{	
			computeBorders_(iter);		
		}

		filterForSize_();
		
		// filter two times for overlaps
		filterForOverlaps_( );
		filterForOverlaps_( );
		
		filterForSignificance_();
					
		cout << iso_map_.size() << " clusters remained after filtering." << endl;
}

void BaseSweepSeeder::voteForCharge_()
{
	// charge states > 10 should rareley be encountered
	vector<ProbabilityType> charge_scores(10, numeric_limits<ProbabilityType>::max() );

	for (TableIteratorType iter = iso_map_.begin(); iter != iso_map_.end(); ++iter)
	{
		charge_scores.clear();
		
		for (std::vector< ScoredChargeType >::const_iterator scmz_iter = iter->second.scored_charges_.begin();
		       scmz_iter != iter->second.scored_charges_.end();
					 ++scmz_iter)
		{
		//cout << "Vote for charge " << scmz_iter->first << " score " << scmz_iter->second << endl;
			
			if ( (scmz_iter->first-1) >= charge_scores.size() || charge_scores.size() == 0)
			{
				charge_scores.resize( scmz_iter->first, numeric_limits<ProbabilityType>::max() );
			}
			
			if ( charge_scores.at( (scmz_iter->first -1) ) > scmz_iter->second )
			{
			charge_scores.at( (scmz_iter->first -1) ) = scmz_iter->second;
			}	
		} // end for ( std::vector< ScoredChargeType > )
		
		//cout << "Done..." << endl;
		
		// search for winning charge
		ProbabilityType max_vote = numeric_limits<ProbabilityType>::max();
		UInt max_charge = 0;
		
		for (UInt i = 0; i < charge_scores.size(); ++i)
		{
			if (charge_scores[i] < max_vote)
			{
				max_charge = (i + 1);
				max_vote     = charge_scores[i];
			}
		
		}
		//cout << "And the winner is " << max_charge << " with score " << max_vote << endl;
		
		iter->second.peaks_.charge_                   = max_charge;
		iter->second.peaks_.max_charge_score_ = max_vote;
	}
}

BaseSweepSeeder::TableIteratorType BaseSweepSeeder::checkInPreviousScans_(const ScoredMZType& sc_mz, const UInt currscan_index)
{
	   // hash entry to write in
    TableIteratorType entry_to_insert;
		CoordinateType mz_in_hash = 0;
		
// 		cout << "checkInPreviousScans_(...) : Retrieving m/z " << sc_mz.first << " " << currscan_index << endl;
		
		CoordinateType curr_mz = traits_->getPeakMz( make_pair(currscan_index,sc_mz.first) );

    if (currscan_index > 0 && iso_map_.size() > 0 )  // Are we in the first scan?
    {
        // there were some isotopic cluster in the last scan...
        TableConstIteratorType table_iter = searchClosestCluster_(curr_mz);
				CoordinateType delta_mz = fabs(table_iter->first - curr_mz);

				#ifdef DEBUG_FEATUREFINDER
				cout << "m/z distance to closest cluster : " << delta_mz << endl;
				#endif
				
        if ( delta_mz > mass_tolerance_cluster_) // check if first peak of last cluster is close enough
        {
            mz_in_hash = curr_mz; // update current hash key

						#ifdef DEBUG_FEATUREFINDER
            cout << "Last peak cluster too far, creating new cluster" << endl;
						cout << "Tolerance : " << mass_tolerance_cluster_ << endl;
						cout << "Creating cluster at m/z " << mz_in_hash << std::endl;
						#endif

            // create new isotopic cluster
            IsotopeClusterScoredCharge isoclust;
            isoclust.scans_.push_back( currscan_index );          
            entry_to_insert = iso_map_.insert( TableType::value_type(mz_in_hash, isoclust) );
        }
        else
        {
						#ifdef DEBUG_FEATUREFINDER
            cout << "Found matching cluster within distance (m/z) " << delta_mz << endl;
						#endif
						
						// there is at least one cluster with m/z within the tolerance, 
						// so have to check how far away (in rt) this cluster is						
            mz_in_hash = table_iter->first;	// set hash key

            pair<TableIteratorType, TableIteratorType> range = iso_map_.equal_range(mz_in_hash);
            
						// we want to find the previous scan
						// currentscan_index can't be zero so we don't have to check for that.
            UInt scan_wanted = (currscan_index - 1);
				
						if (!checkForMatchingCluster_(range,scan_wanted,entry_to_insert))
						{
							// nope, too far
							#ifdef DEBUG_FEATUREFINDER
							cout << "But too far away in rt or in the same scan" << endl;
        			cout << "=> Creating new cluster at m/z: " << curr_mz<< endl;
							#endif
						
							IsotopeClusterScoredCharge isoclust;
        			entry_to_insert = iso_map_.insert( TableType::value_type(curr_mz, isoclust) );						
						}
            
						// save current rt and m/z
            entry_to_insert->second.scans_.push_back( currscan_index );
        }
    }
    else // we are in the first scan (or haven't found any peak cluster), so we have to create a new cluster 
    {
				#ifdef DEBUG_FEATUREFINDER
        cout << "First scan => creating new cluster." << endl;
        cout << "Creating new cluster at m/z: " << curr_mz<< endl;
				#endif

        mz_in_hash = curr_mz; // update current hash key

        IsotopeClusterScoredCharge isoclust;
        isoclust.scans_.push_back( currscan_index );
        entry_to_insert = iso_map_.insert( TableType::value_type(mz_in_hash, isoclust) );

    } // end if (iso_last_scan.size() > 0)

    return entry_to_insert;
}

bool BaseSweepSeeder::checkForMatchingCluster_(const pair<TableIteratorType, TableIteratorType>& range, const UInt scan_wanted, TableIteratorType& entry_to_insert )
{

	// so far we check only for matching cluster in the previous scan
	// maybe it makes sense to take masses being further away into account
	UInt closest_scan = 0;
	bool scan_found   = false;
		
	// loop over all matching cluster
	for (TableIteratorType iter = range.first; iter != range.second; ++iter)
	{
   		// enumerate all scans
    	// the scan number we are searching for is not necessarily the last one
    	// in this cluster if there were other very close local maxima in the same scan.
    	for (vector<UInt>::const_iterator it = iter->second.scans_.begin();
      	     it != iter->second.scans_.end();
            ++it)
    	{
					// find closest scan number
					if (*it >= closest_scan && *it <= scan_wanted)
					{
						closest_scan   = *it;				
						entry_to_insert = iter;		// remember iterator
						scan_found      = true;
					}				
    	}
	} // end for (TableIteratorType )
	
	UInt dist = (scan_wanted - closest_scan);	
	if (dist < rt_tolerance_cluster_ && scan_found)
	{		
		return true;	
	}
	
	// That's too far
	return false;	
}

void BaseSweepSeeder::sumUp_(SpectrumType& scan, UInt current_scan_index)
{	
		for ( UInt i=current_scan_index + 1; i <= current_scan_index + scans_to_sumup_ && i < traits_->getData().size() ; ++i )
    {
				AlignAndSum_(scan,traits_->getData()[i]);
    }
}

void BaseSweepSeeder::substractLastScan_(SpectrumType& current_scan, UInt current_scan_index)
{
		if (current_scan_index >= 1)
		{		
			cout << "Substracting " << (current_scan_index-1) << endl; 
			AlignAndSubstract_(current_scan,traits_->getData()[ current_scan_index-1 ]);		
		}
}


void BaseSweepSeeder::addNextScan_(SpectrumType& scan, UInt current_scan_index)
{
	UInt next_scan_index = current_scan_index + scans_to_sumup_;	

	if (next_scan_index < traits_->getData().size() )
	{
		AlignAndSum_(scan,traits_->getData()[ next_scan_index ]);
	}
}

void BaseSweepSeeder::AlignAndSum_(SpectrumType& scan, const SpectrumType& neighbour)
{
    if (scan.size() == 0 || neighbour.size() == 0)
        return;

    UInt index_newscan = 0;
     for (SpectrumType::const_iterator p = neighbour.begin(); p != neighbour.end(); ++p)
    {
        while (scan[index_newscan].getMZ() < p->getMZ() && index_newscan < scan.size())
            ++index_newscan;

        // This seems to happen more frequently than expected -> quit the loop
        if (index_newscan >= scan.size() )
            break;

        if (index_newscan > 0)
        {
            CoordinateType left_diff   = fabs(scan[index_newscan-1].getMZ() - p->getMZ());
            CoordinateType right_diff = fabs(scan[index_newscan].getMZ() - p->getMZ());

            // check which neighbour is closer
            if (left_diff < right_diff && (left_diff < mass_tolerance_alignment_) )
            {
                scan[ (index_newscan-1) ].setIntensity( scan[ (index_newscan-1) ].getIntensity() + p->getIntensity() );
            }
            else if (right_diff < mass_tolerance_alignment_)
            {
                scan[ (index_newscan) ].setIntensity( scan[ (index_newscan) ].getIntensity() + p->getIntensity() );
            }
        }
        else // no left neighbour available
        {
            CoordinateType right_diff = fabs(scan[index_newscan].getMZ() - p->getMZ());
            if (right_diff < mass_tolerance_alignment_)
            {
                scan[index_newscan].setIntensity( scan[index_newscan].getIntensity() + p->getIntensity() );
            }
        }
    } // end for (all peaks in neighbouring scan)
		
}	// end of AlignAndSum_(....)


void BaseSweepSeeder::AlignAndSubstract_(SpectrumType& scan, const SpectrumType& neighbour)
{
    if (scan.size() == 0 || neighbour.size() == 0)
        return;

    UInt index_newscan = 0;
    for (SpectrumType::const_iterator p = neighbour.begin(); p != neighbour.end(); ++p)
    {

        while (scan[index_newscan].getMZ() < p->getMZ() && index_newscan < scan.size())
					++index_newscan;
        
				// This seems to happen more frequently than expected -> quit the loop
        if (index_newscan >= scan.size() )
            break;

        if (index_newscan > 0)
        {
            CoordinateType left_diff   = fabs(scan[index_newscan-1].getMZ() - p->getMZ());
            CoordinateType right_diff = fabs(scan[index_newscan].getMZ() - p->getMZ());

            // check which neighbour is closer
            if (left_diff < right_diff && (left_diff < mass_tolerance_alignment_) )
            {
                scan[ (index_newscan-1) ].setIntensity( scan[ (index_newscan-1) ].getIntensity() - p->getIntensity() );
            }
            else if (right_diff < mass_tolerance_alignment_)
            {
                scan[ (index_newscan) ].setIntensity( scan[ (index_newscan) ].getIntensity() - p->getIntensity() );
            }
        }
        else // no left neighbour available
        {
            CoordinateType right_diff = fabs(scan[index_newscan].getMZ() - p->getMZ());
            if (right_diff < mass_tolerance_alignment_)
            {
                scan[index_newscan].setIntensity( scan[index_newscan].getIntensity() - p->getIntensity() );
            }
        }
    } // end for (all peaks in neighbouring scan)
		
}	// end of AlignAndSum_(....)



}	// end of namespace OpenMS

