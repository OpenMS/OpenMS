// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_IDMAPPER_H
#define OPENMS_ANALYSIS_ID_IDMAPPER_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <limits>

namespace OpenMS 
{
  /**
    @brief Annotates an MSExperiment, FeatureMap or ConsensusMap with peptide identifications
 		
	  ProteinIdentifications are assigned to the whole map.
	  
 		The retention time and mass-to-charge ratio of the PeptideIdentification have to be given in the MetaInfoInterface as the values "MZ" and "RT".
 		
 		m/z-matching on peptide side can be done either with the precursor m/z value of the peptide identification or with the theoretical masses of the peptide hits (see "mz_reference" parameter).

		See the documentation of the individual @p annotate methods for more in-depth information.
 		
 		@htmlinclude OpenMS_IDMapper.parameters
 		
  */
  class OPENMS_DLLAPI IDMapper
		: public DefaultParamHandler
  {
    public:
			enum Measure {MEASURE_PPM=0, MEASURE_DA};
			
      /// Default constructor
      IDMapper();
      
      /// Copy C'Tor	
      IDMapper(const IDMapper& cp);
      
      /// Assignment
      IDMapper& operator= (const IDMapper& rhs);
			
			/**
				@brief Mapping method for peak maps
				
				The identifications stored in a PeptideIdentification instance can be added to the
    		corresponding spectrum.
				
			  @param map MSExperiment to receive the identifications
			  @param ids PeptideIdentification for the ConsensusFeatures
			  @param protein_ids ProteinIdentification for the ConsensusMap
      	
				@exception Exception::MissingInformation is thrown if the MetaInfoInterface of @p ids does not contain 'MZ' and 'RT'.
      */
      template <typename PeakType>				
      void annotate(MSExperiment<PeakType>& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids)
  		{
				checkHits_(ids);
				
				//append protein identifications
				map.getProteinIdentifications().insert(map.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());
		
  			//store mapping of scan RT to index
				std::multimap<DoubleReal, Size> experiment_precursors;
				for (Size i = 0; i < map.size(); i++)
				{
					experiment_precursors.insert(std::make_pair(map[i].getRT(), i));
				}
				
				//store mapping of identification RT to index
				std::multimap<DoubleReal, Size> identifications_precursors;
				for (Size i = 0; i < ids.size(); i++)
				{
					identifications_precursors.insert(std::make_pair(ids[i].getMetaValue("RT"), i));
				}
				
				//calculate the actual mapping
				std::multimap<DoubleReal, Size>::iterator experiment_iterator = experiment_precursors.begin();
				std::multimap<DoubleReal, Size>::iterator identifications_iterator = identifications_precursors.begin();	
				while(experiment_iterator != experiment_precursors.end() && identifications_iterator != identifications_precursors.end())
				{
					while(identifications_iterator != identifications_precursors.end())
					{
						// testing whether the retention times are within the precision threshold
						if (fabs(experiment_iterator->first - identifications_iterator->first) < rt_tolerance_)
						{
							// testing whether the m/z fits
							if (!map[experiment_iterator->second].getPrecursors().empty())
							{
								if (fabs((DoubleReal)(ids[identifications_iterator->second].getMetaValue("MZ")) - map[experiment_iterator->second].getPrecursors()[0].getMZ()) < mz_tolerance_)
								{
									if (!(ids[identifications_iterator->second].empty()))
									{
										map[experiment_iterator->second].getPeptideIdentifications().push_back(ids[identifications_iterator->second]);
									}
								}
							}
						}
						++identifications_iterator;
					}
					identifications_iterator = identifications_precursors.begin();	
					++experiment_iterator;
				}
  		}
    
			/**
				@brief Mapping method for feature maps

		 		If @em all features have at least one convex hull, peptide positions are matched against the bounding boxes of the convex hulls by default. If not, the positions of the feature centroids are used. The respective coordinates of the centroids are also used for matching (in place of the corresponding ranges from the bounding boxes) if @p use_centroid_rt or @p use_centroid_mz are true.

				In any case, tolerance in RT and m/z dimension is applied according to the global parameters @p rt_tolerance and @p mz_tolerance. Tolerance is understood as "plus or minus x", so the matching range is actually increased by twice the tolerance value.
	
			  If several features (incl. tolerance) overlap the position of a peptide identification, the identification is annotated to all of them.

				@param map FeatureMap to receive the identifications
			  @param ids PeptideIdentification for the ConsensusFeatures
			  @param protein_ids ProteinIdentification for the ConsensusMap
			  @param use_centroid_rt Whether to use the RT value of feature centroids even if convex hulls are present
			  @param use_centroid_mz Whether to use the m/z value of feature centroids even if convex hulls are present
			
				@exception Exception::MissingInformation is thrown if the MetaInfoInterface of @p ids does not contain "MZ" and "RT"
			*/
			template <typename FeatureType>
			void annotate(FeatureMap<FeatureType>& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids, bool use_centroid_rt=false, bool use_centroid_mz=false)
			{
				// std::cout << "Starting annotation..." << std::endl;
				checkHits_(ids);

				// append protein identifications
				map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

				// check if all features have at least one convex hull
				// if not, use the centroid and the given tolerances
				if (!(use_centroid_rt && use_centroid_mz))
				{
					for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
					{
						if (f_it->getConvexHulls().size()==0)
						{
							use_centroid_rt = true;
							use_centroid_mz = true;
							LOG_WARN << "IDMapper warning: at least one feature has no convex hull - using centroid coordinates for matching" << std::endl;
							break;
						}
					}
				}

				bool use_avg_mass = false; // use avg. peptide masses for matching?
				if (use_centroid_mz && (param_.getValue("mz_reference") == "peptide"))
				{
					// if possible, check which m/z value is reported for features,
					// so the appropriate peptide mass can be used for matching
					use_avg_mass = checkMassType_(map.getDataProcessing());
				}

				// calculate feature bounding boxes only once:
				std::vector< DBoundingBox<2> > boxes;				
				DoubleReal min_rt = std::numeric_limits<DoubleReal>::max(), 
					max_rt = -std::numeric_limits<DoubleReal>::max();
				// std::cout << "Precomputing bounding boxes..." << std::endl;
				boxes.reserve(map.size());
				for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); 
						f_it != map.end(); ++f_it)
				{
					DBoundingBox<2> box;
					if (!(use_centroid_rt && use_centroid_mz))
					{
						box = f_it->getConvexHull().getBoundingBox();
					}
					if (use_centroid_rt)
					{
						box.setMinX(f_it->getRT());
						box.setMaxX(f_it->getRT());
					}
					if (use_centroid_mz)
					{
						box.setMinY(f_it->getMZ());
						box.setMaxY(f_it->getMZ());
					}
					increaseBoundingBox_(box);
					boxes.push_back(box);

					min_rt = std::min(min_rt, box.minPosition().getX());
					max_rt = std::max(max_rt, box.maxPosition().getX());
				}

				// hash bounding boxes of features by RT:
				// RT range is partitioned into slices (bins) of 1 second; every feature
				// that overlaps a certain slice is hashed into the corresponding bin
				std::vector<std::vector<SignedSize> > hash_table;
				// make sure the RT hash table has indices >= 0 and doesn't waste space
				// in the beginning:
				SignedSize offset(0);

        if (map.size()>0)
        {
				  // std::cout << "Setting up hash table..." << std::endl;
				  offset = SignedSize(floor(min_rt));
          // this only works if features were found
				  hash_table.resize(SignedSize(floor(max_rt)) - offset + 1);			
				  for (Size index = 0; index < boxes.size(); ++index)
				  {
					  const DBoundingBox<2>& box = boxes[index];
					  for (SignedSize i = SignedSize(floor(box.minPosition().getX())); 
							   i <= SignedSize(floor(box.maxPosition().getX())); ++i)
					  {
						  hash_table[i - offset].push_back(index);
					  }
				  }
        }
        else
        {
          LOG_WARN << "IDMapper received an empty FeatureMap! All peptides are mapped as 'unassigned'!" << std::endl;
        }
		 
				// for statistics:
				Size matches_none = 0, matches_single = 0, matches_multi = 0;
				
				// std::cout << "Finding matches..." << std::endl;
				// iterate over peptide IDs:
				for (std::vector<PeptideIdentification>::const_iterator id_it = 
							 ids.begin(); id_it != ids.end(); ++id_it)
				{
					// std::cout << "Peptide ID: " << id_it - ids.begin() << std::endl;

					if (id_it->getHits().empty()) continue;

					DoubleList mz_values;
					DoubleReal rt_value;
					IntList charges;
					getIDDetails_(*id_it, rt_value, mz_values, charges, use_avg_mass);

					if ((rt_value < min_rt) || (rt_value > max_rt)) // RT out of bounds
					{
 						map.getUnassignedPeptideIdentifications().push_back(*id_it);
						++matches_none;
						continue;
					}

					// iterate over candidate features:
					Size index = SignedSize(floor(rt_value)) - offset;
					Size matching_features = 0;
					for (std::vector<SignedSize>::iterator hash_it = 
								 hash_table[index].begin(); hash_it != hash_table[index].end();
							 ++hash_it)
					{
						Feature& feat = map[*hash_it];

						bool check_charge = use_charge_; // need to check the charge state?
						if (check_charge && (mz_values.size() == 1)) // check now
						{
							if (!charges.contains(feat.getCharge())) continue;
							check_charge = false; // don't need to check later
						}

						// iterate over m/z values (only one if "mz_ref." is "precursor"):
						Size index = 0;
						for (DoubleList::iterator mz_it = mz_values.begin();
								 mz_it != mz_values.end(); ++mz_it, ++index)
						{
							if (check_charge && (charges[index] != feat.getCharge()))
							{
								continue; // charge states need to match
							}

							DPosition<2> id_pos(rt_value, *mz_it);
							if (boxes[*hash_it].encloses(id_pos)) // potential match
							{
								if (use_centroid_mz)
								{ 
									// only one m/z value to check, which was alredy incorporated
									// into the overall bounding box -> success!
									feat.getPeptideIdentifications().push_back(*id_it);
									++matching_features;
									break; // "mz_it" loop
								}
								// else: check all the mass traces
								bool found_match = false;
								for (std::vector<ConvexHull2D>::iterator ch_it = 
											 feat.getConvexHulls().begin(); ch_it != 
											 feat.getConvexHulls().end(); ++ch_it)
								{
									DBoundingBox<2> box = ch_it->getBoundingBox();
									if (use_centroid_rt)
									{
										box.setMinX(feat.getRT());
										box.setMaxX(feat.getRT());
									}
									increaseBoundingBox_(box);
									if (box.encloses(id_pos)) // success!
									{ 
										feat.getPeptideIdentifications().push_back(*id_it);
										++matching_features;
										found_match = true;
										break; // "ch_it" loop
									}
								}
								if (found_match) break; // "mz_it" loop
							}
						}
					}
					if (matching_features == 0)
					{
 						map.getUnassignedPeptideIdentifications().push_back(*id_it);
						++matches_none;
					}
					else if (matching_features == 1) ++matches_single;
					else ++matches_multi;
				}

				//some statistics output
				LOG_INFO << "Unassigned peptides: " << matches_none << "\n"
								 << "Peptides assigned to exactly one feature: " 
								 << matches_single << "\n"
								 << "Peptides assigned to multiple features: " 
								 << matches_multi << std::endl;

			}
			
			/**
				@brief Mapping method for consensus maps

			  If several consensus features lie inside the allowed deviation, the peptide identifications
			  are mapped to all the consensus features.

			  @param map ConsensusMap to receive the identifications
			  @param ids PeptideIdentification for the ConsensusFeatures
			  @param protein_ids ProteinIdentification for the ConsensusMap
			  @param measure_from_subelements Do distance estimate from FeatureHandles instead of Centroid
			 
				@exception Exception::MissingInformation is thrown if the MetaInfoInterface of @p ids does not contain 'MZ' and 'RT'
			*/
		  void annotate(ConsensusMap& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids, bool measure_from_subelements=false);     
		
		protected:
			void updateMembers_();
						
			///Allowed RT deviation
			DoubleReal rt_tolerance_;
			///Allowed m/z deviation
			DoubleReal mz_tolerance_;
			///Measure used for m/z
			Measure measure_;
			///Do charge states have to match?
			bool use_charge_;
			
			/// compute absolute Da tolerance, for a given m/z,
			/// when @p measure is MEASURE_DA, the value is unchanged,
			/// for MEASURE_PPM it is computed according to currently allowed ppm tolerance
			DoubleReal getAbsoluteMZTolerance_(const DoubleReal mz) const;
			
			/// check if distance constraint is fulfilled (using @p rt_tolerance_, @p mz_tolerance_ and @p measure_)
			bool isMatch_(const DoubleReal rt_distance, const DoubleReal mz_theoretical, const DoubleReal mz_observed) const;
			
			/// helper function that checks if all peptide hits are annotated with RT and MZ meta values
			void checkHits_(const std::vector<PeptideIdentification>& ids) const;
			
			/// get RT, m/z and charge value(s) of a PeptideIdentification
			/// - multiple m/z values are returned if "mz_reference" is set to "peptide" (one for each PeptideHit)
			/// - one m/z value is returned if "mz_reference" is set to "precursor"
			void getIDDetails_(const PeptideIdentification& id, DoubleReal& rt_pep, DoubleList& mz_values, IntList& charges, bool use_avg_mass=false) const;

			/// increase a bounding box by the given RT and m/z tolerances
			void increaseBoundingBox_(DBoundingBox<2>& box);

			/// try to determine the type of m/z value reported for features, return
			/// whether average peptide masses should be used for matching
			bool checkMassType_(const std::vector<DataProcessing>& processing) const;

  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDMAPPER_H
