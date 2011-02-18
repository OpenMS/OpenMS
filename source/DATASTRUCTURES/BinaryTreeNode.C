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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>


namespace OpenMS
{

  BinaryTreeNode::BinaryTreeNode(const Size i, const Size j, const Real x) : left_child(i), right_child(j), distance(x)
	{
	}

	BinaryTreeNode::BinaryTreeNode(const BinaryTreeNode& source) : left_child(source.left_child), right_child(source.right_child), distance(source.distance)
	{
	}

	BinaryTreeNode::~BinaryTreeNode()
	{
	}

	BinaryTreeNode& BinaryTreeNode::operator = (const BinaryTreeNode& source)
	{
		if (this != &source)
		{
			left_child = source.left_child;
			right_child = source.right_child;
			distance = source.distance;
		}
		return *this;
	}

}

