// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Steffen Sass $
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/DATASTRUCTURES/DataSubset.h>

namespace OpenMS
{
  DataSubset::DataSubset(DataPoint& data_point)
  {
	  data_points.push_back(&data_point);
	  rt = data_point.rt;
	  mz = data_point.mz;
  }
  DataSubset::DataSubset(const DataSubset& copy) :GridElement(copy)
  {
	  data_points = copy.data_points;
	  tree.insert(tree.begin(),copy.tree.begin(),copy.tree.end());
	  distance_iterators = copy.distance_iterators;
  }
  DataSubset::DataSubset(const DataSubset* copy_ptr)
  {
	  data_points=copy_ptr->data_points;
	  tree.insert(tree.begin(),copy_ptr->tree.begin(),copy_ptr->tree.end());
	  distance_iterators=copy_ptr->distance_iterators;
	  mz=copy_ptr->mz;
	  rt=copy_ptr->rt;
  }
  Int DataSubset::operator<(const DataSubset &el) const
  {
	  std::list<DataPoint*> data1=this->data_points;
	  std::list<DataPoint*> data2=el.data_points;
	  return data1.size() < data2.size();
  }

  bool DataSubset::operator !=(const DataSubset &el) const
  {
	  if (this->data_points.size()!=el.data_points.size())
	  {
		  return true;
	  }
	  else
	  {
		  std::list<DataPoint*> data1=this->data_points;
		  std::list<DataPoint*>::iterator it1=data1.begin();
		  std::list<DataPoint*> data2=el.data_points;
		  std::list<DataPoint*>::iterator it2=data2.begin();
		  while((*it1)->mz==(*it2)->mz && (*it1)->rt==(*it2)->rt)
		  {
			  ++it1;
			  ++it2;
			  if (it1==data1.end() && it2==data2.end())
				  return false;
		  }
		  return true;
	  }
  }

  bool DataSubset::operator ==(const DataSubset & el) const
  {
	  return !(*this!=el);
  }

  int DataSubset::size()
  {
	  return data_points.size();
  }

  int DataSubset::getID() const
  {
	  return data_points.front()->feature_id;
  }
}
