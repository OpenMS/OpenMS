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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut, Steffen Sass $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/SILACTreeNode.h>
#include <cmath>

namespace OpenMS
{
  SILACTreeNode::SILACTreeNode()
  {

  }

  SILACTreeNode::SILACTreeNode(DataPoint* data1_, DataPoint* data2_, DoubleReal distance_)
  {
    data1 = data1_;
    data2 = data2_;
    distance = distance_;
  }

  bool SILACTreeNode::operator == (const SILACTreeNode &cp) const
  {
	  if( this->data1 != cp.data1) return false;
	  if( this->data2 != cp.data2) return false;
	  if( this->distance != cp.distance) return false;
	  return true;
  }

  bool SILACTreeNode::operator != (const SILACTreeNode &cp) const
  {
    return !(*this == cp);
  }

  bool SILACTreeNode::operator < (const SILACTreeNode &cp) const
  {
	  if (std::abs(this->distance - cp.distance) <= 0.00000001)
    {
		  return *(this->data1) < *(cp.data1);
    }
	  else
    {
		  return this->distance < cp.distance;
    }
  }
}

