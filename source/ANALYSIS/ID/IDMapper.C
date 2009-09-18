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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMapper.h>

using namespace std;

namespace OpenMS 
{
  IDMapper::IDMapper()
  	: rt_delta_(5.0),
  		mz_delta_(0.01)
  { 
  }

			
	DoubleReal IDMapper::getRTDelta() const
	{
		return rt_delta_;
	}
	
	void IDMapper::setRTDelta(DoubleReal rt_delta)
	{
		rt_delta_ = rt_delta;
	}
	
	DoubleReal IDMapper::getMZDelta() const
	{
		return mz_delta_;
	}
	
	void IDMapper::setMZDelta(DoubleReal mz_delta)
	{
		mz_delta_ = mz_delta;
	}

  void IDMapper::annotate(ConsensusMap& map, const std::vector<PeptideIdentification>& ids, const std::vector<ProteinIdentification>& protein_ids, bool measure_from_subelements)
	{
		checkHits_(ids);
				
		//append protein identifications to Map
		map.getProteinIdentifications().insert(map.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());

		//keep track of assigned/unassigned peptide identifications
		std::set<Size> assigned;
					
		// store which peptides fit which feature (and avoid double entries)
		// consensusMap -> {peptide_index}
		std::vector < std::set< size_t> > mapping(map.size());
		
		//iterate over the peptide IDs
		for (Size i=0; i<ids.size(); ++i)
		{
			if (ids[i].getHits().size()==0) continue;

			DoubleReal rt_pep = ids[i].getMetaValue("RT");
			DoubleReal mz_pep = ids[i].getMetaValue("MZ");

			//iterate over the features
			for(Size cm_index = 0 ; cm_index<map.size(); ++cm_index)
			{
				//check if we compare distance from centroid or subelements
				if (!measure_from_subelements)
				{
					if ( (fabs(rt_pep-map[cm_index].getRT()) <= rt_delta_) && (fabs(mz_pep-map[cm_index].getMZ()) <= mz_delta_)  )
					{
						map[cm_index].getPeptideIdentifications().push_back(ids[i]);
						assigned.insert(i);
					}
				}
				else
				{
					for(ConsensusFeature::HandleSetType::const_iterator it_handle = map[cm_index].getFeatures().begin(); 
							it_handle != map[cm_index].getFeatures().end(); 
							++it_handle)
					{
						if ( (fabs(rt_pep - it_handle->getRT()) <= rt_delta_) && (fabs(mz_pep - it_handle->getMZ()) <= mz_delta_) )
						{
							if (mapping[cm_index].count(i) == 0)
							{
								map[cm_index].getPeptideIdentifications().push_back(ids[i]);
								assigned.insert(i);
								mapping[cm_index].insert(i);
							}
							continue; // we added this peptide already.. no need to check further
						}
					}
				}
			}
		}

		//append unassigned peptide identifications
		for (Size i=0; i<ids.size(); ++i)
		{
			if (assigned.count(i)==0)
			{
				map.getUnassignedPeptideIdentifications().push_back(ids[i]);
			}
		}

	}

} // namespace OpenMS
