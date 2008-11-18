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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDConsensusFeatureMapper.h>

using namespace std;

namespace OpenMS 
{
  IDConsensusFeatureMapper::IDConsensusFeatureMapper() 
  {
  }
  
  void IDConsensusFeatureMapper::annotate(ConsensusMap& map,
																					const std::vector<PeptideIdentification>& ids,
																					const std::vector<ProteinIdentification>& protein_ids,
																					DoubleReal mz_delta,
																					DoubleReal rt_delta,
																					bool measure_from_subelements)
	{		
		//append protein identifications to Map
		map.getProteinIdentifications().insert(map.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());

		
		// store which peptides fit which feature (and avoid double entries)
		// consensusMap -> {peptide_index}
		std::vector < std::set< size_t> > mapping(map.size());
		
		//iterate over the peptide IDs
		for (size_t i=0; i<ids.size(); ++i)
		{
			if (ids[i].getHits().size()==0)
			{
				continue;
			}

			if (!ids[i].metaValueExists("RT"))
			{
				throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDConsensusFeatureMapper: meta data value 'RT' missing for peptide identification!"); 
			}

			if (!ids[i].metaValueExists("MZ"))
			{
				throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDConsensusFeatureMapper: meta data value 'MZ' missing for peptide identification!"); 
			}

			DoubleReal rt_pep = ids[i].getMetaValue("RT");
			DoubleReal mz_pep = ids[i].getMetaValue("MZ");

			//iterate over the features
			for(size_t cm_index = 0 ; cm_index<map.size(); ++cm_index)
			{
				//check if we compare distance from centroid or subelements
				if (!measure_from_subelements)
				{
					if ( (fabs(rt_pep-map[cm_index].getRT()) <= rt_delta) && (fabs(mz_pep-map[cm_index].getMZ()) <= mz_delta)  )
					{
						map[cm_index].getPeptideIdentifications().push_back(ids[i]);
					}
				}
				else
				{
					for(ConsensusFeature::HandleSetType::iterator it_handle = map[cm_index].getFeatures().begin(); 
							it_handle != map[cm_index].getFeatures().end(); 
							++it_handle)
					{
						if ( (fabs(rt_pep - it_handle->getRT()) <= rt_delta) && (fabs(mz_pep - it_handle->getMZ()) <= mz_delta) )
						{
							if (mapping[cm_index].count(i) == 0)
							{
								map[cm_index].getPeptideIdentifications().push_back(ids[i]);
								mapping[cm_index].insert(i);
							}
							continue; // we added this peptide already.. no need to check further
						}
					}
				}

			} // ! features loop
		} // ! peptides loop

	}
} // namespace OpenMS
