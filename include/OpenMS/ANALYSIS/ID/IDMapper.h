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
      IDMapper& operator = (const IDMapper& rhs);
			
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
				checkHits_(ids);
				
				//append protein identifications
				map.getProteinIdentifications().insert(map.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());

				//check if all feature have at least one convex hull
				//if not, use the centroid and the given deltas
				if (!use_centroids)
				{
					for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
					{
						if (f_it->getConvexHulls().size()==0)
						{
							use_centroids = true;
							std::cout << "IDMapper warning: at least one feature has no convex hull => using centroids!" << std::endl;
							break;
						}
					}
				}
				
				//precalculate feature bounding boxes
				std::vector< DBoundingBox<2> > bbs;
				if (!use_centroids)
				{
					bbs.reserve(map.size());
					for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
					{
						DBoundingBox<2> bb = f_it->getConvexHull().getBoundingBox();
						bb.setMin(bb.min() - DPosition<2>(rt_delta_, getAbsoluteMZDelta_(bb.min().getY())));
						bb.setMax(bb.max() + DPosition<2>(rt_delta_, getAbsoluteMZDelta_(bb.max().getY())));
						bbs.push_back(bb);
					}
				}
				
				//keep track of assigned/unassigned peptide identifications
				std::map<Size, Size> assigned;
				
				std::vector<DoubleReal> mz_values;
				DoubleReal rt_pep;
			
				// features...
				std::vector< DBoundingBox<2> >::const_iterator bb_it = bbs.begin();
				for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
				{
					//iterate over the IDs
					for (Size i=0; i<ids.size(); ++i)
					{
						if (ids[i].getHits().size()==0) continue;

						getRTandMZofID_(ids[i], rt_pep, mz_values);
						// if set to TRUE, we leave the i_mz-loop as we added the whole ID with all hits
						bool was_added=false; // was current pep-m/z matched?!

						// iterate over m/z values of pepIds
						for (Size i_mz=0;i_mz<mz_values.size();++i_mz)
						{
							DoubleReal mz_pep = mz_values[i_mz];

							if (use_centroids)
							{
							

								if ( isMatch_(rt_pep - f_it->getRT(), mz_pep, f_it->getMZ()) )
								{
									was_added = true;
									f_it->getPeptideIdentifications().push_back(ids[i]);
									assigned[i]++;
								}
							}
							else
							{
								DPosition<2> id_pos(rt_pep, mz_pep);
								//check if the ID lies within the bounding box

								if (bb_it->encloses(id_pos))
								{
									// iterate over all convex hulls
									for(std::vector<ConvexHull2D>::iterator ch_it = f_it->getConvexHulls().begin(); ch_it!=f_it->getConvexHulls().end(); ++ch_it)
									{
										DBoundingBox<2> bb = ch_it->getBoundingBox();
										bb.setMin(bb.min() - DPosition<2>(rt_delta_, getAbsoluteMZDelta_(bb.min().getY())));
										bb.setMax(bb.max() + DPosition<2>(rt_delta_, getAbsoluteMZDelta_(bb.max().getY())));
										if (bb.encloses(id_pos))
										{
											was_added = true;
											f_it->getPeptideIdentifications().push_back(ids[i]);
											assigned[i]++;
											break;
										}
									}
								} 
							} // !centroids
							
							if (was_added) break;
						} // m/z values to check

					} // ID's
					if(!use_centroids)	++bb_it;
				} // features
				
				Size matches_none = 0;
				Size matches_single = 0;
				Size matches_multi = 0;
				
				//append unassigned peptide identifications
				for (Size i=0; i<ids.size(); ++i)
				{
					Size matches = assigned[i];
					if (matches==0)
					{
						map.getUnassignedPeptideIdentifications().push_back(ids[i]);
						matches_none++;
					}
					else if (matches==1)
					{
						matches_single++;
					}
					else
					{
						matches_multi++;
					}
				}
				
				//some statistics output
				std::cout << "Unassigned peptides: " << matches_none << std::endl;
				std::cout << "Peptides assigned to exactly one features: " << matches_single << std::endl;
				std::cout << "Peptides assigned to multiple features: " << matches_multi << std::endl;
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
			
			///get RT and M/Z value(s) of a peptideIdentification
			/// - multiple m/z values are returned if "mz_reference" is set to "PeptideMass" (one for each PeptideHit)
			/// - one m/z value is returned if "mz_reference" is set to "PrecursorMZ"
			void getRTandMZofID_(const PeptideIdentification& id, DoubleReal& rt_pep, std::vector<DoubleReal>& mz_values) const;
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDMAPPER_H
