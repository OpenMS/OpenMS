// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_BINARYTREENODE_H
#define OPENMS_DATASTRUCTURES_BINARYTREENODE_H

#include <OpenMS/CONCEPT/Types.h>


namespace OpenMS
{
  
	/** @brief Elements of a binary tree used to represent a hierarchical clustering process

			strict indexing/topology is assumed, i.e. node no. x represents clusteringstep no. x
			left_child and right_child are each the lowest indices to elements of the merged clusters, distance is the distance of the two children
	*/
	class OPENMS_DLLAPI BinaryTreeNode
	{
		public:
		/// constructor
		BinaryTreeNode(const Size i, const Size j, const Real x);

		/// destructor
		~BinaryTreeNode();

		/// copy constructor
		BinaryTreeNode(const BinaryTreeNode& source);

		/// assignment operator
		BinaryTreeNode& operator = (const BinaryTreeNode& source);

		Size left_child;
		Size right_child;
		Real distance;

		private:
		BinaryTreeNode();
	};

}

#endif /* BINARYTREENODE_H_ */
