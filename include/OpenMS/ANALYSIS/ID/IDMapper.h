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
// $Maintainer: $
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
    @brief Annotates a MSExperiment, FeatureMap or ConsensusMap with peptide identifications
 		
	  ProteinIdentifications are assigned to the whole map.
	  
 		The retention time and mass-to-charge ratio of the PeptideIdentification have to 
 		be given in the MetaInfoInterface as the values 'MZ' and 'RT'.
  */
  class OPENMS_DLLAPI IDMapper
  {
    public:

      /// Default constructor
      IDMapper();
			
			/// Returns the allowed RT deviation (default: 5.0)
			DoubleReal getRTDelta() const;
			/// Sets the allowed RT deviation
			void setRTDelta(DoubleReal rt_delta);
			
			/// Returns the allowed RT deviation (default: 0.01)
			DoubleReal getMZDelta() const;
			/// Sets the allowed RT deviation
			void setMZDelta(DoubleReal mz_delta);

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

		 		If @em all features have at least one convex hull, the identification are mapped to the contex hull.
		 		If not, the allowed m/z and RT deviation from the feature centroid position is checked.
		
			  If several features lie inside the allowed deviation, the peptide identifications
			  are mapped to all the features.

				@param map FeatureMap to receive the identifications
			  @param ids PeptideIdentification for the ConsensusFeatures
			  @param protein_ids ProteinIdentification for the ConsensusMap
			  @param use_centroids If @em true, the feature centoids are used, even if convex hulls are present
			
				@exception Exception::MissingInformation is thrown if the MetaInfoInterface of @p ids does not contain 'MZ' and 'RT'
			*/
			template <typename FeatureType>
		  void annotate(FeatureMap<FeatureType>& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids, bool use_centroids=false)
			{
				checkHits_(ids);
				
				//append protein identifications
				map.getProteinIdentifications().insert(map.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());

				//check if all feature have at least one convex hull
				//if not, use the controid and the given deltas
				if (!use_centroids)
				{
					for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
					{
						if (f_it->getConvexHulls().size()==0)
						{
							use_centroids = true;
							std::cout << "IDMapper warning: at lease one feature has no convex hull => using centroids!" << std::endl;
							break;
						}
					}
				}
				
				//precalculate feature bounding boxes
				DPosition<2> delta (rt_delta_, mz_delta_);
				std::vector< DBoundingBox<2> > bbs;
				if (!use_centroids)
				{
					bbs.reserve(map.size());
					for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
					{
						DBoundingBox<2> bb = f_it->getConvexHull().getBoundingBox();
						bb.setMin(bb.min() - delta);
						bb.setMax(bb.max() + delta);
						bbs.push_back(bb);
					}
				}
				
				//keep track of assigned/unassigned peptide identifications
				std::map<Size, Size> assigned;
				std::vector< DBoundingBox<2> >::const_iterator bb_it = bbs.begin();
				for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
				{
					//iterate over the IDs
					for (Size i=0; i<ids.size(); ++i)
					{
						if (ids[i].getHits().size()==0) continue;

						if (use_centroids)
						{
							if ( (fabs((DoubleReal)ids[i].getMetaValue("RT")-f_it->getRT()) <= rt_delta_) 
								&& (fabs((DoubleReal)ids[i].getMetaValue("MZ")-f_it->getMZ()) <= mz_delta_)  )
							{
								f_it->getPeptideIdentifications().push_back(ids[i]);
								assigned[i]++;
							}
						}
						else
						{
							DPosition<2> id_pos(ids[i].getMetaValue("RT"),ids[i].getMetaValue("MZ"));
							//check if the ID lies within the bounding box

							if (bb_it->encloses(id_pos))
							{
								// iterate over all convex hulls
								for(std::vector<ConvexHull2D>::iterator ch_it = f_it->getConvexHulls().begin(); ch_it!=f_it->getConvexHulls().end(); ++ch_it)
								{
									DBoundingBox<2> bb = ch_it->getBoundingBox();
									bb.setMin(bb.min() - delta);
									bb.setMax(bb.max() + delta);
									if (bb.encloses(id_pos))
									{
										f_it->getPeptideIdentifications().push_back(ids[i]);
										assigned[i]++;
										break;
									}
								}
							}
						}
					}
					if(!use_centroids)	++bb_it;
				}
				
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
			
			///Allowed RT deviation
			DoubleReal rt_delta_;
			///Allowed m/z deviation
			DoubleReal mz_delta_;
			
			///Helper function that checks if all peptide hits are annotated with RT and MZ meta values
			void checkHits_(const std::vector<PeptideIdentification>& ids)
			{
				for (Size i=0; i<ids.size(); ++i)
				{
					if (!ids[i].metaValueExists("RT"))
					{
						throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDMapper: meta data value 'RT' missing for peptide identification!"); 
					}
					if (!ids[i].metaValueExists("MZ"))
					{
						throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDMapper: meta data value 'MZ' missing for peptide identification!"); 
					}
				}
			}
			
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDMAPPER_H
