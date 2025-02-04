// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ML/CLUSTERING/AverageLinkage.h>
#include <OpenMS/ML/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <vector>
//#include <iostream>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AverageLinkage, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AverageLinkage* ptr = nullptr;
AverageLinkage* nullPointer = nullptr;
START_SECTION(AverageLinkage())
{
	ptr = new AverageLinkage();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~AverageLinkage())
{
	delete ptr;
}
END_SECTION

START_SECTION((AverageLinkage(const AverageLinkage &source)))
{
	AverageLinkage al1;
	AverageLinkage copy(al1);
}
END_SECTION

START_SECTION((AverageLinkage& operator=(const AverageLinkage &source)))
{
	AverageLinkage copy,al2;
	copy = al2;
}
END_SECTION

START_SECTION((void operator()(DistanceMatrix< float > &original_distance, std::vector<BinaryTreeNode>& cluster_tree, const float threshold=1) const))
{
	DistanceMatrix<float> matrix(6,666);
	matrix.setValue(1,0,0.5f);
	matrix.setValue(2,0,0.8f);
	matrix.setValue(2,1,0.3f);
	matrix.setValue(3,0,0.6f);
	matrix.setValue(3,1,0.8f);
	matrix.setValue(3,2,0.8f);
	matrix.setValue(4,0,0.8f);
	matrix.setValue(4,1,0.8f);
	matrix.setValue(4,2,0.8f);
	matrix.setValue(4,3,0.4f);
	matrix.setValue(5,0,0.7000001f); //~ minimal adjustment for gcc 4 with -o2
	matrix.setValue(5,1,0.8f);
	matrix.setValue(5,2,0.8f);
	matrix.setValue(5,3,0.8f);
	matrix.setValue(5,4,0.8f);
	DistanceMatrix<float> matrix2(matrix);

	vector< BinaryTreeNode > result;
	vector< BinaryTreeNode > tree;
	//~ tree.push_back(BinaryTreeNode(1,2,0.3f));
	//~ tree.push_back(BinaryTreeNode(2,3,0.4f));
	//~ tree.push_back(BinaryTreeNode(0,1,0.65f));
	//~ tree.push_back(BinaryTreeNode(0,1,0.766667f));
	//~ tree.push_back(BinaryTreeNode(0,1,0.78f));
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(3,4,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.65f));
	tree.push_back(BinaryTreeNode(0,3,0.766667f));
	tree.push_back(BinaryTreeNode(0,5,0.78f));

	AverageLinkage al;
	al(matrix,result);
	TEST_EQUAL(tree.size(), result.size());
	for (Size i = 0; i < result.size(); ++i)
	{
			TEST_EQUAL(tree[i].left_child, result[i].left_child);
			TEST_EQUAL(tree[i].right_child, result[i].right_child);
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(tree[i].distance, result[i].distance);
	}

	float th(0.7f);
	tree.pop_back();
	tree.pop_back();
	tree.push_back(BinaryTreeNode(0,3,-1.0f));
	tree.push_back(BinaryTreeNode(0,5,-1.0f));
	result.clear();

	al(matrix2,result,th);
	TEST_EQUAL(tree.size(), result.size());
	for (Size i = 0; i < result.size(); ++i)
	{
			TEST_EQUAL(tree[i].left_child, result[i].left_child);
			TEST_EQUAL(tree[i].right_child, result[i].right_child);
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_REAL_SIMILAR(tree[i].distance, result[i].distance);
	}

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



