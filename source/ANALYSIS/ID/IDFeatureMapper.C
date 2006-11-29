// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/KERNEL/DimensionDescription.h>

using namespace std;

namespace OpenMS 
{
  IDFeatureMapper::IDFeatureMapper() 
  {
    
  }
  
  void IDFeatureMapper::annotate(DFeatureMap<2> fm, const vector<Identification>& identifications, const vector<float>& precursor_retention_times, const vector<float>& precursor_mz_values)
  	throw (Exception::Precondition)
	{
		//Precondition
		if (precursor_retention_times.size()!=precursor_mz_values.size() || precursor_retention_times.size()!=identifications.size())
		{
			throw Exception::Precondition(__FILE__,__LINE__,__PRETTY_FUNCTION__,"Identification, RT and m/z vectos must have the same size!");
		}
		
		//Dimensions
    enum DimensionId
    {
      RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
      MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };
		
		//iterate over the features
		for(DFeatureMap<2>::Iterator f_it = fm.begin(); f_it!=fm.end(); ++f_it)
		{
			cout << endl << "* Feature (rt/mz): " << f_it->getPosition()[RT] << " " << f_it->getPosition()[MZ] << endl;
			DBoundingBox<2> bb = f_it->getBoundingBox();
			const DFeature<2>::ConvexHullVector& ch_vec = f_it->getConvexHulls();
			
			//iterate over the IDs
			for (UnsignedInt i=0; i<precursor_retention_times.size(); ++i)
			{
				cout << "  * ID (rt/mz): " << precursor_retention_times[i] << " " << precursor_mz_values[i] << endl;
				DPosition<2> id_pos(precursor_retention_times[i],precursor_mz_values[i]);
				
				//check if the ID lies within the bouning box. if it does not => next id
				if (!bb.encloses(id_pos)) continue;
				cout << "  * ID inside BB " << endl;
				
				for(DFeature<2>::ConvexHullVector::const_iterator ch_it = ch_vec.begin(); ch_it!=ch_vec.end(); ++ch_it)
				{
					cout << "    * Convex Hull" << endl;
					if (ch_it->encloses(id_pos))
					{
						cout << "    * !!HIT!!" << endl;
						continue;
					}
				}		
			}

		}
	}
} // namespace OpenMS
