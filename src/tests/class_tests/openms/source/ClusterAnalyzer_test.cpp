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
#include <OpenMS/ML/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ClusterAnalyzer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ClusterAnalyzer* ptr = nullptr;
ClusterAnalyzer* nullPointer = nullptr;
START_SECTION(ClusterAnalyzer())
{
	ptr = new ClusterAnalyzer();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ClusterAnalyzer())
{
	delete ptr;
}
END_SECTION

START_SECTION((ClusterAnalyzer(const ClusterAnalyzer &source)))
{
  NOT_TESTABLE
}
END_SECTION

ptr = new ClusterAnalyzer();

START_SECTION((std::vector< float > averageSilhouetteWidth(const std::vector< BinaryTreeNode > &tree, const DistanceMatrix< float > &original)))
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
	matrix.setValue(5,0,0.7f);
	matrix.setValue(5,1,0.8f);
	matrix.setValue(5,2,0.8f);
	matrix.setValue(5,3,0.8f);
	matrix.setValue(5,4,0.8f);

	vector<float> asw(5);
	asw[0]=0.170833f;
	asw[1]=0.309722f;
	asw[2]=0.306412f;
	asw[3]=0.125744f;
	asw[4]=0;

	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(3,4,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,3,0.6f));
	tree.push_back(BinaryTreeNode(0,5,0.7f));
	vector<float> result = ptr->averageSilhouetteWidth(tree, matrix);
	TEST_EQUAL(result.size(), asw.size());
	for (Size i = 0; i < result.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.001);
			TEST_REAL_SIMILAR(result[i], asw[i]);
	}

}
END_SECTION

START_SECTION((std::vector< float > dunnIndices(const std::vector<BinaryTreeNode>& tree, const DistanceMatrix<float>& original, const bool tree_from_singlelinkage = false)))
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
	matrix.setValue(5,0,0.7f);
	matrix.setValue(5,1,0.8f);
	matrix.setValue(5,2,0.8f);
	matrix.setValue(5,3,0.8f);
	matrix.setValue(5,4,0.8f);

  vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(3,4,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,3,0.6f));
	tree.push_back(BinaryTreeNode(0,5,0.7f));
	vector<float> di(5);
	di[0]=0.4f/0.3f; di[1]=0.5f/0.4f; di[2]=0.6f/0.8f; di[3]=0.7f/0.8f; di[4]=0.0f;
	vector<float> result = ptr->dunnIndices(tree, matrix);
	TEST_EQUAL(result.size(), di.size());
	for (Size i = 0; i < result.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.001);
			TEST_REAL_SIMILAR(result[i], di[i]);
	}
	result = ptr->dunnIndices(tree, matrix, true);
	TEST_EQUAL(result.size(), di.size());
	for (Size i = 0; i < result.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.001);
			TEST_REAL_SIMILAR(result[i], di[i]);
	}
}
END_SECTION

START_SECTION((std::vector< float > cohesion(const std::vector< std::vector<Size> >& clusters, const DistanceMatrix<float>& original)))
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
	matrix.setValue(5,0,0.7f);
	matrix.setValue(5,1,0.8f);
	matrix.setValue(5,2,0.8f);
	matrix.setValue(5,3,0.8f);
	matrix.setValue(5,4,0.8f);

  Size a[] = {0,1,2,3,4,5};
	vector< vector<Size> > clusters;
	clusters.push_back(vector<Size>(a,a+3));
	clusters.push_back(vector<Size>(a+3,a+5));
	clusters.push_back(vector<Size>(a+5,a+6));

	vector<float> cohesions;
	cohesions.push_back(0.533f);
	cohesions.push_back(0.4f);
	cohesions.push_back(0.7f);

	vector<float> result = ptr->cohesion(clusters, matrix);
	TEST_EQUAL(cohesions.size(), result.size());
	for (Size i = 0; i < cohesions.size(); ++i)
	{
		TOLERANCE_ABSOLUTE(0.001);
		TEST_REAL_SIMILAR(cohesions[i], result[i]);
	}

	clusters.clear();
	clusters.push_back(vector<Size>(a,a+4));
	clusters.push_back(vector<Size>(a+4,a+5));
	clusters.push_back(vector<Size>(a+5,a+6));

	cohesions.clear();
	cohesions.push_back(0.633f);
	cohesions.push_back(0.7f);
	cohesions.push_back(0.7f);

	result = ptr->cohesion(clusters, matrix);
	TEST_EQUAL(cohesions.size(), result.size());
	for (Size i = 0; i < cohesions.size(); ++i)
	{
		TOLERANCE_ABSOLUTE(0.001);
		TEST_REAL_SIMILAR(cohesions[i], result[i]);
	}
}
END_SECTION

START_SECTION((float averagePopulationAberration(Size cluster_quantity, std::vector<BinaryTreeNode>& tree)))
{
	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(3,4,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,3,0.6f));
	tree.push_back(BinaryTreeNode(0,5,0.7f));

	float result = ptr->averagePopulationAberration(3, tree);
	TEST_REAL_SIMILAR(2.0/3.0, result);
}
END_SECTION

START_SECTION((void cut(const Size cluster_quantity, const std::vector<BinaryTreeNode>& tree, std::vector<std::vector<Size> >& clusters)))
{
  Size a[] = {0,1,2,3,4,5};
	vector< vector<Size> > clusters;
	vector< vector<Size> > result;

	result.push_back(vector<Size>(a,a+3));
	result.push_back(vector<Size>(a+3,a+5));
	result.push_back(vector<Size>(a+5,a+6));

  vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(3,4,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,3,0.6f));
	tree.push_back(BinaryTreeNode(0,5,0.7f));
  ptr->cut(3, tree, clusters);
	TEST_EQUAL(clusters.size(), result.size());
	for (Size i = 0; i < clusters.size(); ++i)
	{
		TEST_EQUAL(clusters[i].size(), result[i].size());
		for (Size j = 0; j < clusters[i].size(); ++j)
		{
			TEST_EQUAL(clusters[i][j], result[i][j]);
		}
	}

	Size b[] = {0,1,5,8,10,12,2,3,9,11,4,6,7};
	result.clear();

	result.push_back(vector<Size>(b,b+1));
	result.push_back(vector<Size>(b+1,b+6));
	result.push_back(vector<Size>(b+6,b+10));
	result.push_back(vector<Size>(b+10,b+13));

  std::vector< BinaryTreeNode > trunk;
	trunk.push_back(BinaryTreeNode(4,6,0.1f));
	trunk.push_back(BinaryTreeNode(2,3,0.11f));
	trunk.push_back(BinaryTreeNode(5,8,0.111f));
	trunk.push_back(BinaryTreeNode(4,7,0.2f));
	trunk.push_back(BinaryTreeNode(2,9,0.22f));
	trunk.push_back(BinaryTreeNode(1,10,0.222f));
	trunk.push_back(BinaryTreeNode(2,11,0.3f));
	trunk.push_back(BinaryTreeNode(1,5,0.33f));
	trunk.push_back(BinaryTreeNode(1,12,0.333f));
	trunk.push_back(BinaryTreeNode(0,1,-1.0f));
	trunk.push_back(BinaryTreeNode(0,2,-1.0f));
	trunk.push_back(BinaryTreeNode(0,4,-1.0f));
	clusters.clear();
  ptr->cut(4, trunk, clusters);
	TEST_EQUAL(clusters.size(), result.size());
	for (Size i = 0; i < clusters.size(); ++i)
	{
		TEST_EQUAL(clusters[i].size(), result[i].size());
		for (Size j = 0; j < clusters[i].size(); ++j)
		{
			TEST_EQUAL(clusters[i][j], result[i][j]);
		}
	}

}
END_SECTION

START_SECTION((void cut(const Size cluster_quantity, const std::vector<BinaryTreeNode>& tree, std::vector< std::vector<BinaryTreeNode> >& subtrees)))
{
	std::vector< std::vector< BinaryTreeNode > > c_ts(4),ts;
  std::vector< BinaryTreeNode > trunk;
	trunk.push_back(BinaryTreeNode(4,6,0.1f));
	trunk.push_back(BinaryTreeNode(2,3,0.11f));
	trunk.push_back(BinaryTreeNode(5,8,0.111f));
	trunk.push_back(BinaryTreeNode(4,7,0.2f));
	trunk.push_back(BinaryTreeNode(2,9,0.22f));
	trunk.push_back(BinaryTreeNode(1,10,0.222f));
	trunk.push_back(BinaryTreeNode(2,11,0.3f));
	trunk.push_back(BinaryTreeNode(1,5,0.33f));
	trunk.push_back(BinaryTreeNode(1,12,0.333f));
	trunk.push_back(BinaryTreeNode(0,1,-1.0f));
	trunk.push_back(BinaryTreeNode(0,2,-1.0f));
	trunk.push_back(BinaryTreeNode(0,4,-1.0f));

	//~ c_ts[0].push_back(BinaryTreeNode(0,1,-1.0f));
	//~ c_ts[0].push_back(BinaryTreeNode(0,2,-1.0f));
	//~ c_ts[0].push_back(BinaryTreeNode(0,4,-1.0f));
	c_ts[1].push_back(BinaryTreeNode(5,8,0.111f));
	c_ts[1].push_back(BinaryTreeNode(1,10,0.222f));
	c_ts[1].push_back(BinaryTreeNode(1,5,0.33f));
	c_ts[1].push_back(BinaryTreeNode(1,12,0.333f));
	c_ts[2].push_back(BinaryTreeNode(2,3,0.11f));
	c_ts[2].push_back(BinaryTreeNode(2,9,0.22f));
	c_ts[2].push_back(BinaryTreeNode(2,11,0.3f));
	c_ts[3].push_back(BinaryTreeNode(4,6,0.1f));
	c_ts[3].push_back(BinaryTreeNode(4,7,0.2f));

	ptr->cut(4, trunk, ts);
	TEST_EQUAL(ts.size(), c_ts.size());
	for (Size i = 0; i < c_ts.size() && i < ts.size(); ++i)
	{
		TEST_EQUAL(ts[i].size(), c_ts[i].size());
		for (Size j = 0;  j < ts[i].size() && j < c_ts[i].size(); ++j)
		{
			TEST_EQUAL(ts[i][j].right_child, c_ts[i][j].right_child);
			TEST_EQUAL(ts[i][j].left_child, c_ts[i][j].left_child);
			TEST_EQUAL(ts[i][j].distance, c_ts[i][j].distance);
		}
	}
}
END_SECTION

START_SECTION((String newickTree(const std::vector<BinaryTreeNode>& tree, const bool include_distance = false)))
{
	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(3,4,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,3,0.6f));
	tree.push_back(BinaryTreeNode(0,5,0.7f));

	String result = ptr->newickTree(tree);
	TEST_EQUAL(result,"( ( ( 0 , ( 1 , 2 ) ) , ( 3 , 4 ) ) , 5 )");
	result = ptr->newickTree(tree,true);
	TEST_EQUAL(result,"( ( ( 0:0.5 , ( 1:0.3 , 2:0.3 ):0.5 ):0.6 , ( 3:0.4 , 4:0.4 ):0.6 ):0.7 , 5:0.7 )");
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



