// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm $
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
	  
 		The retention time and mass-to-charge ratio of the PeptideIdentification have to 
 		be given in the MetaInfoInterface as the values 'MZ' and 'RT'.
 		
 		m/z-Matching on the peptide side can be done either with the precursor m/z value of the peptideIdentification or
 		the masses of the peptideHits (see 'mz_reference' parameter).
 		
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
						if (fabs(experiment_iterator->first - identifications_iterator->first) < rt_delta_)
						{
							// testing whether the m/z fits
							if (!map[experiment_iterator->second].getPrecursors().empty())
							{
								if (fabs((DoubleReal)(ids[identifications_iterator->second].getMetaValue("MZ")) - map[experiment_iterator->second].getPrecursors()[0].getMZ()) < mz_delta_)
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

		 		If @em all features have at least one convex hull, the identifications are mapped to the convex hull.
		 		If not, the allowed m/z and RT deviation from the feature centroid (RT,MZ) position is checked.
		
			  If several features lie inside the allowed deviation, the peptide identifications
			  are mapped to all the features.

				@param map FeatureMap to receive the identifications
			  @param ids PeptideIdentification for the ConsensusFeatures
			  @param protein_ids ProteinIdentification for the ConsensusMap
			  @param use_centroids If @em true, the feature centroids (RT,MZ position) are used, even if convex hulls are present
			
				@exception Exception::MissingInformation is thrown if the MetaInfoInterface of @p ids does not contain 'MZ' and 'RT'
			*/
			template <typename FeatureType>
		  void annotate(FeatureMap<FeatureType>& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids, bool use_centroids=false)
			{
				// std::cout << "Starting annotation..." << std::endl;
				checkHits_(ids);
				
				// append protein identifications
				map.getProteinIdentifications().insert(map.getProteinIdentifications().end(), protein_ids.begin(), protein_ids.end());

				// check if all features have at least one convex hull
				// if not, use the centroid and the given deltas
				if (!use_centroids)
				{
					for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
					{
						if (f_it->getConvexHulls().size()==0)
						{
							use_centroids = true;
							LOG_WARN << "IDMapper warning: at least one feature has no convex hull => using centroids!" << std::endl;
							break;
						}
					}
				}
				
				// hash features (bounding boxes) by RT:
				// RT range is partitioned into slices (bins) of 1 second; every feature
				// that overlaps a certain slice is hashed into the corresponding bin
				// std::cout << "Setting up hash table..." << std::endl;
				std::vector<std::vector<SignedSize> > hash_table;
				DoubleReal min_rt = std::numeric_limits<DoubleReal>::max(), 
					max_rt = -std::numeric_limits<DoubleReal>::max();
				// make sure the RT hash table has indices > 0	
				SignedSize offset;

				// calculate feature bounding boxes only once (if applicable):
				std::vector< DBoundingBox<2> > boxes;
				if (use_centroids)
				{
					// fill the hash table:
					map.sortByRT();
					min_rt = map.front().getRT() - rt_delta_;
					max_rt = map.back().getRT() + rt_delta_;
					offset = SignedSize(floor(min_rt));
					hash_table.resize(SignedSize(floor(max_rt)) - offset + 1);			
					for (Size index = 0; index < map.size(); ++index)
					{
						DoubleReal feat_rt = map[index].getRT();
						for (SignedSize i = SignedSize(floor(feat_rt - rt_delta_)); 
								 i <= SignedSize(floor(feat_rt + rt_delta_)); ++i)
						{
							hash_table[i - offset].push_back(index);
						}
					}
				}
				else // use bounding boxes
				{
					// std::cout << "Precomputing bounding boxes..." << std::endl;
					boxes.reserve(map.size());
					for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); 
							f_it != map.end(); ++f_it)
					{
						DBoundingBox<2> box;
						if (param_.getValue("mz_reference") == "PeptideMass")
						{ 
							if (f_it->getConvexHulls().size() == 1)
							{ // only one hull for the whole feature (-> "isotope_wavelet"):
								box = f_it->getConvexHull().getBoundingBox();
								box.setMinY(f_it->getMZ());
								box.setMaxY(f_it->getMZ());
							}
							else
							{ // find monoisotopic mass trace:
								std::vector<ConvexHull2D>::iterator mono_it = 
									min_element(f_it->getConvexHulls().begin(), 
															f_it->getConvexHulls().end(),
															IDMapper::mass_trace_comp_);
								box = mono_it->getBoundingBox();
							}
						}
						else
						{
							box = f_it->getConvexHull().getBoundingBox();
						}
						increaseBoundingBox_(box);
						boxes.push_back(box);

						// "min"/"max" redefinition (DBoundingBox) confuses the compiler:
						// min_rt = min(min_rt, box.min().getX());
						// max_rt = max(max_rt, box.max().getX());
						if (box.min().getX() < min_rt) min_rt = box.min().getX();
						if (box.max().getX() > max_rt) max_rt = box.max().getX();
					}
				
					// fill the hash table:
					// std::cout << "Filling hash table..." << std::endl;
					offset = SignedSize(floor(min_rt));
					hash_table.resize(SignedSize(floor(max_rt)) - offset + 1);			
					for (Size index = 0; index < boxes.size(); ++index)
					{
						const DBoundingBox<2>& box = boxes[index];
						for (SignedSize i = SignedSize(floor(box.min().getX())); 
								 i <= SignedSize(floor(box.max().getX())); ++i)
						{
							hash_table[i - offset].push_back(index);
						}
					}
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

					std::vector<DoubleReal> mz_values;
					DoubleReal rt_value;
					getRTandMZofID_(*id_it, rt_value, mz_values);

					if ((rt_value < min_rt) || (rt_value > max_rt)) // RT out of bounds
					{
 						map.getUnassignedPeptideIdentifications().push_back(*id_it);
						++matches_none;
						continue;
					}

					// iterate over candidate features:
					Size index = SignedSize(floor(rt_value)) - offset;
					Size matching_features = 0;
					for (std::vector<SignedSize>::iterator hash_it = hash_table[index].begin();
							 hash_it != hash_table[index].end(); ++hash_it)
					{
						Feature& feat = map[*hash_it];

						// iterate over m/z values
						// (only one if "mz_reference" is "PrecursorMZ"):
						for (std::vector<DoubleReal>::iterator mz_it = mz_values.begin();
								 mz_it != mz_values.end(); ++mz_it)
						{
							if (use_centroids)
							{
								if (isMatch_(rt_value - feat.getRT(), *mz_it, feat.getMZ()))
								{
									feat.getPeptideIdentifications().push_back(*id_it);
									++matching_features;
									break; // "mz_it" loop
								}
							}
							else // use bounding boxes
							{
								DPosition<2> id_pos(rt_value, *mz_it);
								if (boxes[*hash_it].encloses(id_pos))
								{
									if (param_.getValue("mz_reference") == "PeptideMass")
									{ // already checked monoisotopic mass trace -> success!
										feat.getPeptideIdentifications().push_back(*id_it);
										++matching_features;
										break; // "mz_it" loop
									}
									// else: check all the mass traces
									// (in this case, "mz_values" contains only one value ->
									// no need to break out of the "mz_it" loop)
									for(std::vector<ConvexHull2D>::iterator ch_it = 
												feat.getConvexHulls().begin(); ch_it != 
												feat.getConvexHulls().end(); ++ch_it)
									{
										DBoundingBox<2> box = ch_it->getBoundingBox();
										increaseBoundingBox_(box);
										if (box.encloses(id_pos))
										{ // success!
											feat.getPeptideIdentifications().push_back(*id_it);
											++matching_features;
											break; // "ch_it" loop
										}
									}
								}
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
								 << "Peptides assigned to exactly one feature: "  << matches_single << "\n"
								 << "Peptides assigned to multiple features: "  << matches_multi 
								 << std::endl;
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
			DoubleReal rt_delta_;
			///Allowed m/z deviation
			DoubleReal mz_delta_;
			///Measure used for m/z
			Measure measure_;
			
			/// compute absolute Da delta, for a given m/z,
			/// when @p measure is MEASURE_DA, the value is unchanged,
			/// for MEASURE_PPM it is computed according to currently allowed ppm delta
			DoubleReal getAbsoluteMZDelta_(const DoubleReal mz) const;
			
			/// check if distance constraint is fulfilled (using @p rt_delta_, @p mz_delta_ and @p measure_)
			bool isMatch_(const DoubleReal rt_distance, const DoubleReal mz_theoretical, const DoubleReal mz_observed) const;
			
			///Helper function that checks if all peptide hits are annotated with RT and MZ meta values
			void checkHits_(const std::vector<PeptideIdentification>& ids) const;
			
			///get RT and M/Z value(s) of a PeptideIdentification
			/// - multiple m/z values are returned if "mz_reference" is set to "PeptideMass" (one for each PeptideHit)
			/// - one m/z value is returned if "mz_reference" is set to "PrecursorMZ"
			void getRTandMZofID_(const PeptideIdentification& id, DoubleReal& rt_pep, std::vector<DoubleReal>& mz_values) const;

			/// "operator<" to compare mass traces (convex hulls) by mean m/z
			static bool mass_trace_comp_(const ConvexHull2D& first, 
																	 const ConvexHull2D& second);

			/// increase a bounding box by the given RT and m/z deltas
			void increaseBoundingBox_(DBoundingBox<2>& box);

  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDMAPPER_H
