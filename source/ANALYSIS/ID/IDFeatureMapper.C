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

#include <OpenMS/ANALYSIS/ID/IDFeatureMapper.h>

using namespace std;

namespace OpenMS 
{
  IDFeatureMapper::IDFeatureMapper() 
  {
    
  }
  
  void IDFeatureMapper::annotate(FeatureMap<>& fm, const vector<PeptideIdentification>& ids, const vector<ProteinIdentification>& protein_ids)
	{		
		//append protein identifications
		fm.getProteinIdentifications().insert(fm.getProteinIdentifications().end(),protein_ids.begin(),protein_ids.end());
		
		//iterate over the features
		for(FeatureMap<>::Iterator f_it = fm.begin(); f_it!=fm.end(); ++f_it)
		{
#ifdef DEBUG_ID_MAPPING
			cout << endl << "IDFeatureMapper -- Feature " << f_it->getRT() << " " << f_it->getMZ() << endl;
#endif
			DBoundingBox<2> bb = f_it->getConvexHull().getBoundingBox();
			vector<ConvexHull2D>& ch_vec = f_it->getConvexHulls();
			
			//iterate over the IDs
			for (UInt i=0; i<ids.size(); ++i)
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

#ifdef DEBUG_ID_MAPPING
				cout << "  * ID (rt/mz): " << ids[i].getMetaValue("RT") << " " << ids[i].getMetaValue("MZ") << endl;
#endif
				DPosition<2> id_pos(ids[i].getMetaValue("RT"),ids[i].getMetaValue("MZ"));
				
				//check if the ID lies within the bouning box. if it does not => next id
				if (!bb.encloses(id_pos))
				{
#ifdef DEBUG_ID_MAPPING
					cout << "  * outside BB " << endl;
#endif
					continue;
				}
				for(vector<ConvexHull2D>::iterator ch_it = ch_vec.begin(); ch_it!=ch_vec.end(); ++ch_it)
				{
					if (ch_it->encloses(id_pos))
					{
						f_it->getPeptideIdentifications().push_back(ids[i]);
#ifdef DEBUG_ID_MAPPING
						cout << "    * !!HIT!!" << endl;
#endif
						break;
					}
				}		
			}
		}
	}
} // namespace OpenMS
