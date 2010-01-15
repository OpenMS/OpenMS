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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConsensusMap.h>

using namespace std;

namespace OpenMS
{

#if 0
	bool ConsensusMap::isValid(String& error_message) const
	{
		//store the size of the input files
		Map<Size,Size> map_sizes;
		for (FileDescriptions::const_iterator it=file_description_.begin(); it!=file_description_.end(); ++it)
    {
    	if (it->second.size==0)
    	{
    		map_sizes[it->first] = numeric_limits<Size>::max();
			}
			else
			{
    		map_sizes[it->first] = it->second.size;
			}
		}

		//actual check
		for (Size i = 0; i < size(); ++i)
    {
      const ConsensusFeature& f = operator[](i);
      for (ConsensusFeature::HandleSetType::const_iterator it = f.begin(); it!=f.end(); ++it)
      {
      	//check if all map indices are registered
      	if (!file_description_.has(it->getMapIndex()))
      	{
      		error_message = String("ConsensusFeature ") + i + ": Invalid map index " + it->getMapIndex();
      		return false;
      	}
      	//check if the element indices are valid
      	if (it->getElementIndex()>=map_sizes[it->getMapIndex()])
      	{
      		error_message = String("ConsensusFeature ") + i + ": Invalid element index " + it->getElementIndex() + " (size is " + map_sizes[it->getMapIndex()] + ")";
      		return false;
      	}
      }
    }
    return true;
	}
#endif

  ostream& operator << (ostream& os, const ConsensusMap& cons_map)
  {
		for (ConsensusMap::FileDescriptions::const_iterator it=cons_map.getFileDescriptions().begin(); it!=cons_map.getFileDescriptions().end(); ++it)
    {
    	os << "Map " << it->first << ": " << it->second.filename << " - " << it->second.label << " - " << it->second.size << endl;
    }

    for (Size i = 0; i < cons_map.size(); ++i)
    {
      os << cons_map[i] << endl;
    }

    return os;
  }

	void ConsensusMap::updateRanges()
	{
		clearRanges();
		updateRanges_(begin(),end());

		//enlarge the range by the internal points of each feature
		for (Size i=0; i<size(); ++i)
		{
			for (ConsensusFeature::HandleSetType::const_iterator it=operator[](i).begin(); it!=operator[](i).end(); ++it)
			{
				DoubleReal rt = it->getRT();
				DoubleReal mz = it->getMZ();
				DoubleReal intensity = it->getIntensity();

				//update RT
				if (rt < pos_range_.min()[Peak2D::RT])
				{
					pos_range_.setMinX(rt);
				}
				if (rt > pos_range_.max()[Peak2D::RT])
				{
					pos_range_.setMaxX(rt);
				}
				//update m/z
				if (mz < pos_range_.min()[Peak2D::MZ])
				{
					pos_range_.setMinY(mz);
				}
				if (mz > pos_range_.max()[Peak2D::MZ])
				{
					pos_range_.setMaxY(mz);
				}
				//update intensity
				if (intensity <  int_range_.minX())
				{
					int_range_.setMinX(intensity);
				}
				if (intensity > int_range_.maxX())
				{
					int_range_.setMaxX(intensity);
				}
			}
		}
	}

}
