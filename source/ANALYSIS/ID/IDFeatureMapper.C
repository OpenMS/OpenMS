// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
  
  void IDFeatureMapper::annotate(FeatureMap<>& fm, const vector<IdentificationData>& ids, const vector<ProteinIdentification>& protein_ids)
  	throw (Exception::Precondition)
	{		
		//assign protein identifications
		fm.setProteinIdentifications(protein_ids);
		
		//iterate over the features
		for(FeatureMap<>::Iterator f_it = fm.begin(); f_it!=fm.end(); ++f_it)
		{
			//cout << endl << "* Feature (rt/mz): " << f_it->getRT() << " " << f_it->getMZ() << endl;
			DBoundingBox<2> bb = f_it->getBoundingBox();
			const Feature::ConvexHullVector& ch_vec = f_it->getConvexHulls();
			
			//iterate over the IDs
			for (UnsignedInt i=0; i<ids.size(); ++i)
			{
				//cout << "  * ID (rt/mz): " << ids[i].rt << " " << ids[i].mz << endl;
				DPosition<2> id_pos(ids[i].rt,ids[i].mz);
				
				//check if the ID lies within the bouning box. if it does not => next id
				if (!bb.encloses(id_pos))
				{
					//cout << "  * outside BB " << endl;
					continue;
				}
				for(Feature::ConvexHullVector::const_iterator ch_it = ch_vec.begin(); ch_it!=ch_vec.end(); ++ch_it)
				{
					//cout << "    * Convex Hull" << endl;
					if (ch_it->encloses(id_pos))
					{
						f_it->getIdentifications().push_back(ids[i].id);
						//cout << "    * !!HIT!!" << endl;
						break;
					}
				}		
			}
		}
	}
} // namespace OpenMS
