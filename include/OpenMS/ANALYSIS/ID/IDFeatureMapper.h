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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_IDFEATUREMAPPER_H
#define OPENMS_ANALYSIS_ID_IDFEATUREMAPPER_H

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <vector>

namespace OpenMS 
{
  /**
    @brief Annotates a FeatureMap with peptide identifications
 		
 		If @em all features have at least one convex hull, the identification are mapped to the contex hull.
 		If not, the allowed m/z and RT deviation from the feature centroid position is checked.

	  If several consensus features lie inside the allowed deviation, the peptide identifications
	  are mapped to all the features.
 		
	  ProteinIdentifications are assigned to the whole map.
	  
 		The retention time and mass-to-charge ratio of the PeptideIdentification have to 
 		be given in the MetaInfoInterface ('MZ' and 'RT').
  */
  class IDFeatureMapper
  {
    public:

      /// Default constructor
      IDFeatureMapper();
      
			/**
				@brief This method does the actual mapping
			
				@param map FeatureMap to receive the identifications
			  @param ids PeptideIdentification for the ConsensusFeatures
			  @param protein_ids ProteinIdentification for the ConsensusMap
			  @param mz_delta Allowed m/z deviation from feature centroid position
			  @param rt_delta Allowed RT deviation from feature centroid position
			  @param use_delta If @em true, the m/z and RT deviation parameters are used even if convex hulls are present
			
				@exception Exception::MissingInformation is thrown if the MetaInfoInterface of @p ids does not contain 'MZ' and 'RT'
			*/
			template <typename FeatureType>
		  void annotate(FeatureMap<FeatureType>& map,
		                const std::vector<PeptideIdentification>& ids,
		                const std::vector<ProteinIdentification>& protein_ids,
                    DoubleReal rt_delta = 0.5,
                    DoubleReal mz_delta = 0.05,
                    bool use_delta = false)
			{				
				//append protein identifications
				map.getProteinIdentifications().insert(map.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());
				
				//check if all feature have at least one convex hull
				//if not, use the controid and the given deltas
				if (!use_delta)
				{
					for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
					{
						if (f_it->getConvexHulls().size()==0)
						{
							use_delta = true;
							break;
						}
					}
				}
				
				for(typename FeatureMap<FeatureType>::Iterator f_it = map.begin(); f_it!=map.end(); ++f_it)
				{
					//iterate over the IDs
					for (UInt i=0; i<ids.size(); ++i)
					{
						if (ids[i].getHits().size()==0)
						{
							continue;
						}
						if (!ids[i].metaValueExists("RT"))
						{
							throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDFeatureMapper: meta data value 'RT' missing for peptide identification!"); 
						}
						if (!ids[i].metaValueExists("MZ"))
						{
							throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDFeatureMapper: meta data value 'MZ' missing for peptide identification!"); 
						}
						
						if (use_delta)
						{
							if ( (fabs((DoubleReal)ids[i].getMetaValue("RT")-f_it->getRT()) <= rt_delta) 
								&& (fabs((DoubleReal)ids[i].getMetaValue("MZ")-f_it->getMZ()) <= mz_delta)  )
							{
								f_it->getPeptideIdentifications().push_back(ids[i]);
							}
						}
						else
						{
							DPosition<2> id_pos(ids[i].getMetaValue("RT"),ids[i].getMetaValue("MZ"));
							
							//check if the ID lies within the bounding box
							if (f_it->getConvexHull().getBoundingBox().encloses(id_pos))
							{
								// iterate over all convex hulls
								for(std::vector<ConvexHull2D>::iterator ch_it = f_it->getConvexHulls().begin(); ch_it!=f_it->getConvexHulls().end(); ++ch_it)
								{
									if (ch_it->encloses(id_pos))
									{
										f_it->getPeptideIdentifications().push_back(ids[i]);
										break;
									}
								}
							}
						}
					}
				}
			}
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDFEATUREMAPPER_H
