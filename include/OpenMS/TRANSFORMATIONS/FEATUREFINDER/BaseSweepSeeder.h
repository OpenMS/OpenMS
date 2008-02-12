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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASESWEEPSEEDER_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASESWEEPSEEDER_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BaseSeeder.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>

#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>

#include <map>
#include <limits>

namespace OpenMS
{
  /** 
  	@brief Abstract base class for sweepline seeding modules.
		
		This is a base class for all seeding modules based on (or implementing) the sweepline paradigm
		for peptide quantification / feature detection in LC/MS maps.
		
		All classes derived from this base class perform the following steps (they differ in the way they find isotopic patterns only): Several scans are combined to improve signal-to-noise level. Several signals that are found in close proximity are then combined to one signal. This is done for all spectra. A signal is only accepted when it is found in several consecutive spectra.
		
		@note Scores for charge estimates should be >= 0 where a low score indicates a high confidence e.g. score should be some sort of p value.
		@note Method updateMembers() in each base class should call updateMembers() in this class before doing anything else. 
		
		</table>
		 <tr><td></td><td></td><td>min_number_scans</td>
		 <td>lower bound for the number of scans in which a isotopic pattern must occur 
		  before it is accepted as seeding region</td></tr>
			 <tr><td></td><td></td><td>max_number_scans</td>
		 <td>upper bound for the number of scans in which a isotopic pattern must occur 
		  (used to remove solvents)</td></tr>
		 <tr><td></td><td></td><td>min_number_peaks</td>
		 <td>min. number of data points for a seeding region</td></tr>
		 <tr><td></td><td></td><td>scans_to_sumup</td>
		 <td>number of scans used during alignment </td></tr>
		 <tr><td></td><td></td><td>mass_tolerance_alignment</td>
		 <td>tolerance in m/z for aligned points </td></tr>
		 	<tr><td></td><td></td><td>mass_tolerance_cluster</td>
		 <td>tolerance in m/z during assembly of isotopic point cluster </td></tr>
		 	<tr><td></td><td></td><td>rt_tolerance_cluster</td>
		 <td>tolerance in m/z during assembly of isotopic point cluster</td></tr>
		</table>		
   
		@improvement Seeder should pass a list of good charge states to model fitter instead of only the best one.
		
    @ingroup FeatureFinder
  */
  class BaseSweepSeeder
  	: public BaseSeeder
  {

	  public:	  	
		
			/// intensity of a peak
			typedef FeaFiTraits::IntensityType IntensityType;
			/// coordinate ( in rt or m/z )
		  typedef FeaFiTraits::CoordinateType CoordinateType;
			/// score
		  typedef DoubleReal ProbabilityType;	
			
			/// charge state estimate with associated score
			typedef std::pair<UInt, ProbabilityType> ScoredChargeType;
			/// m/z position in spectrum with charge estimate and score
			typedef std::pair< UInt, ScoredChargeType > ScoredMZType;
			/// container of scored m/z positions
			typedef std::vector< ScoredMZType > ScoredMZVector;
							
			/// Stores information about isotopic cluster plus scored charge estimates in each scan
  		struct IsotopeClusterScoredCharge : IsotopeCluster
  		{
  			  	
    		IsotopeClusterScoredCharge()
      		:	IsotopeCluster(),
						scored_charges_(),
						first_scan_(),
						last_scan_()
				{	
				}
				
    		/// vector of scored charges
    		std::vector< ScoredChargeType > scored_charges_;
				/// first scan
				UInt first_scan_;
				/// last scan
				UInt last_scan_;
			};		

			/// LC-MS map
			typedef FeaFiTraits::MapType MapType;
			/// a single MS spectrum
			typedef MapType::SpectrumType SpectrumType;
			/// a peak
			typedef	 MapType::PeakType PeakType;
			/// hash table mapping m/z values to groups of (isotopic) peaks
			typedef std::multimap<CoordinateType, IsotopeClusterScoredCharge> TableType;	
			/// Iterator in hash table
			typedef TableType::iterator TableIteratorType;
			/// Const iterator in hash table
			typedef TableType::const_iterator TableConstIteratorType;
									
			
		   /// default constructor 
	    BaseSweepSeeder();
	
	    /// copy constructor 
	    BaseSweepSeeder(const BaseSweepSeeder& source);
	
	    /// destructor 
	    virtual ~BaseSweepSeeder();
	
	    /// assignment operator 
	    virtual BaseSweepSeeder& operator = (const BaseSweepSeeder& source);
		
	    /// return next seed 
		 virtual ChargedIndexSet nextSeed() throw (NoSuccessor);
						 		 			   
		protected:
		
			virtual void updateMembers_();		
			
			/// detects isotopic pattern in a scan
			virtual ScoredMZVector detectIsotopicPattern_(SpectrumType& scan ) = 0;
						
			/// Sweeps through scans and detects isotopic patterns
		  void sweep_();
							
			/// Sums the intensities in adjacent scans
		  void sumUp_(SpectrumType& scan, UInt current_scan_index);
			
			/// Aligns two scans and increases intensities of peaks in @p scan if those peaks are present in @p neighbour
			void AlignAndSum_(SpectrumType& scan, const SpectrumType& neighbour);
			
			/// Aligns two scans and substracts the intensities of matching points.
			void AlignAndSubstract_(SpectrumType& scan, const SpectrumType& neighbour);
			
			/// Substracts the last scan
			void substractLastScan_(SpectrumType& scan, UInt current_scan_index);
						
			/// Aligns current scan with its successor
			void addNextScan_(SpectrumType& scan, UInt current_scan_index);
			
			/// Filter hash of point cluster
			void filterHash_();
			
			/// Decide about most likely charge state by majority vote
			void voteForCharge_();
			
			/// Finds the neighbour of the peak denoted by @p current_mz in the previous scan and returns the mass distance
	    TableConstIteratorType searchClosestCluster_(const CoordinateType current_mz)
	    {
	      // perform binary search to find the closest peak cluster in m/z
	      // 	lower_bound finds the first element whose key is not less than current_mz.first
	      TableConstIteratorType iter = iso_map_.lower_bound( current_mz );
	
	      // the peak found by lower_bound does not have to be the closest one, therefore we have
	      // to check both neighbours
	      if ( iter == iso_map_.end() ) // we are at the end and have only one choice
	      {
	      	--iter;
	      }
        // if the found peak is at the beginning of the spectrum,
        // there is not much we can do.
        else if ( iter != iso_map_.begin() )
        {
					// check which neighbour (left or right) is closer
          if ( (iter->first - current_mz) < (current_mz - (--iter)->first) )
          {
          	++iter;    // peak to the right is closer
          }
	      }
				return iter;
	    }
			
			/// check for cluster in previous scans
			TableIteratorType checkInPreviousScans_(const ScoredMZType&,  const UInt);
			
			/// check for matching cluster. called when there are several point cluster with similar masses.
			bool checkForMatchingCluster_(const std::pair<TableIteratorType, TableIteratorType>&, const UInt, TableIteratorType&);
			
			/// computes the median scan number for a hash entry
			void computeBorders_(TableIteratorType& entry);
			
			/// filters sweepline hash for overlapping point cluster with the same charge
			void filterForOverlaps_();
			
			/// filters sweepline hash for tiny (and probably insignificant) regions
			void filterForSize_();
			
			/// filters sweepline hash for regions with low p-value
			void filterForSignificance_();
			
			/// Deletes the hash entries in @p entries.
			void deleteHashEntries_(std::vector<TableIteratorType>& entries);
							
			/// Maps m/z to sets of peaks 
		  TableType iso_map_;
						
			/// Pointer to the current region
			TableType::const_iterator curr_region_;
		
		  /// Indicates whether the extender has been initialized
			bool is_initialized_;
				
			/// Mass tolerance during scan alignment
			CoordinateType mass_tolerance_alignment_;
			
			/// Number of scans used during alignment
			UInt scans_to_sumup_;
			
			/// Mass tolerance for sweepline
			CoordinateType mass_tolerance_cluster_;
			
			/// Rt tolerance for sweepline		
			UInt rt_tolerance_cluster_;
			
			/// Max. distance in rt for merged peak cluster
			UInt max_rt_dist_merging_;
			
			/// Max. distance in mz for merged peak cluster
			CoordinateType max_mz_dist_merging_;
		
  };
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_BASESWEEPSEEDER_H
