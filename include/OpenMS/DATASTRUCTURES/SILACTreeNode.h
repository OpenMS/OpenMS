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


#ifndef OPENMS_DATASTRUCTURES_SILACTREENODE_H
#define OPENMS_DATASTRUCTURES_SILACTREENODE_H

#include <OpenMS/DATASTRUCTURES/DataPoint.h>


namespace OpenMS
{

/**
		@brief A node of an hierarchical clustering tree.
		@ingroup Datastructures
	*/

class OPENMS_DLLAPI SILACTreeNode 
{
  public:
    /** @brief the first data point of the node

	  */
	  DataPoint* data1;

    /** @brief the second data point of the node

	  */
	  DataPoint* data2;

    /** @brief distance between the two points

	  */
	  DoubleReal distance;

    /** @brief default constructor

    */
	  SILACTreeNode();

    /** @brief detailed constructor

			  @param data1_ first data point
			  @param data2_ second data point
			  @param distance_ distance between the data points
	  */
	  SILACTreeNode(DataPoint* data1_,DataPoint* data2_,DoubleReal distance_);
	  bool operator==(const SILACTreeNode &cp) const;
	  bool operator!=(const SILACTreeNode &cp) const;
	  bool operator<(const SILACTreeNode &cp) const;

};

} // namespace

#endif /* SILACTreeNode_H_ */
