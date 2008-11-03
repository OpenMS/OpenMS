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
  
  void IDConsensusFeatureMapper::annotate(ConsensusMap& cm,
																					const std::vector<PeptideIdentification>& ids,
																					const std::vector<ProteinIdentification>& protein_ids,
																					CoordinateType mz_delta,
																					CoordinateType rt_delta,
																					bool measure_from_subelements)
	{		
		//append protein identifications to Map
		cm.getProteinIdentifications().insert(cm.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());

		
		// store which peptides fit which feature (and avoid double entries)
		// consensusMap -> {peptide_index}
		std::vector < std::set< size_t> > mapping(cm.size());
		
		//iterate over the peptide IDs
		for (size_t i=0; i<ids.size(); ++i)
		{
			if (ids[i].getHits().size()==0)
			{
				continue;
			}

			if (!ids[i].metaValueExists("RT"))
			{
				throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDSpectrumMapper: MetaValue 'RT' missing!"); 
			}

			if (!ids[i].metaValueExists("MZ"))
			{
				throw Exception::MissingInformation(__FILE__,__LINE__,__PRETTY_FUNCTION__, "IDSpectrumMapper: MetaValue 'MZ' missing!"); 
			}

			CoordinateType rt_pep = ids[i].getMetaValue("RT");
			CoordinateType mz_pep = ids[i].getMetaValue("MZ");

			#ifdef DEBUG_MAPPER
			std::cout << "PEP rt:" << rt_pep << " mz:" << mz_pep << "\n";
			#endif
			
			//iterate over the features
			for(size_t cm_index = 0 ; cm_index<cm.size(); ++cm_index)
			{
				#ifdef DEBUG_MAPPER
				std::cout << "CF rt:" << cm[cm_index].getRT() << " mz:" << cm[cm_index].getMZ() << "\n";
				#endif
				
				//check if we compare distance from centroid or subelements
				if (!measure_from_subelements)
				{
					if ( (fabs(rt_pep-cm[cm_index].getRT()) <= rt_delta) &&
							 (fabs(mz_pep-cm[cm_index].getMZ()) <= mz_delta)  )
					{
						#ifdef DEBUG_MAPPER
						std::cout << "!!! HIT\n";
						#endif
						if (mapping[cm_index].count(i) == 0)
						{
							cm[cm_index].getPeptideIdentifications().push_back(ids[i]);
						}
						mapping[cm_index].insert(i); //add peptide #i to feature #cm_index
					}
				}
				else
				{
					for(ConsensusFeature::HandleSetType::iterator it_handle = cm[cm_index].getFeatures().begin(); 
							it_handle != cm[cm_index].getFeatures().end(); 
							++it_handle)
					{
						if ( (fabs(rt_pep - it_handle->getRT()) <= rt_delta) &&
							   (fabs(mz_pep - it_handle->getMZ()) <= mz_delta) )
						{
							#ifdef DEBUG_MAPPER
							std::cout << "!!! HIT\n";
							#endif
							if (mapping[cm_index].count(i) == 0)
							{
								cm[cm_index].getPeptideIdentifications().push_back(ids[i]);
							}
							mapping[cm_index].insert(i); //add peptide #i to feature #cm_index
							continue; // we added this peptide already.. no need to check further
						}
					}
				}

			} // ! features loop
		} // ! peptides loop

	}
} // namespace OpenMS
