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
// $Maintainer: Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SingleLinkage, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SingleLinkage* ptr = 0;
START_SECTION(SingleLinkage())
{
	ptr = new SingleLinkage();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION(~SingleLinkage())
{
	delete ptr;
}
END_SECTION

ptr = new SingleLinkage();

START_SECTION((SingleLinkage(const SingleLinkage &source)))
{
	SingleLinkage copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
}
END_SECTION

START_SECTION((SingleLinkage& operator=(const SingleLinkage &source)))
{
	SingleLinkage copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
}
END_SECTION

START_SECTION((void cluster(const DistanceMatrix< Real > &original_distance, std::vector<BinaryTreeNode>& cluster_tree, const Real threshold=1) const))
{

	DistanceMatrix<Real> matrix(6,666);
	matrix.setValue(1,0,0.5);
	matrix.setValue(2,0,0.8);
	matrix.setValue(2,1,0.3);
	matrix.setValue(3,0,0.6);
	matrix.setValue(3,1,0.8);
	matrix.setValue(3,2,0.8);
	matrix.setValue(4,0,0.8);
	matrix.setValue(4,1,0.8);
	matrix.setValue(4,2,0.8);
	matrix.setValue(4,3,0.4);
	matrix.setValue(5,0,0.7);
	matrix.setValue(5,1,0.8);
	matrix.setValue(5,2,0.8);
	matrix.setValue(5,3,0.8);
	matrix.setValue(5,4,0.8);

	vector< BinaryTreeNode > result;
	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3));
	tree.push_back(BinaryTreeNode(2,3,0.4));
	tree.push_back(BinaryTreeNode(0,1,0.5));
	tree.push_back(BinaryTreeNode(0,1,0.6));
	tree.push_back(BinaryTreeNode(0,1,0.7));

	(*ptr).cluster(matrix,result);
	TEST_EQUAL(tree.size(), result.size());
	for (UInt i = 0; i < tree.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_EQUAL(tree[i].left_child, result[i].left_child);
			TEST_EQUAL(tree[i].right_child, result[i].right_child);
			TEST_REAL_SIMILAR(tree[i].distance, result[i].distance);
	}

	//~ Real th(0.6);
	//~ tree.pop_back();
	//~ tree.pop_back();
	//~ result.clear();

	//~ (*ptr).cluster(matrix,result,th);
	//~ for (UInt i = 0; i < tree.size(); ++i)
	//~ {
			//~ TOLERANCE_ABSOLUTE(0.0001);
			//~ TEST_EQUAL(tree[i].left_child, result[i].left_child);
			//~ TEST_EQUAL(tree[i].right_child, result[i].right_child);
			//~ TEST_REAL_SIMILAR(tree[i].distance, result[i].distance);
	//~ }
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(ptr->getProductName(), "SingleLinkage")
}
END_SECTION

START_SECTION((static ClusterFunctor* create()))
{
	ClusterFunctor* cf = SingleLinkage::create();
	SingleLinkage sl;
	TEST_EQUAL(cf->getName(), sl.getName())
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



