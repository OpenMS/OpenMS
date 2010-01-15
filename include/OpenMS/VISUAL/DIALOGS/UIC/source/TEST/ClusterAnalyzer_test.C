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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ClusterAnalyzer, "$Id: ClusterAnalyzer_test.C 4999 2009-03-31 07:05:10Z herrzweiblum $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ClusterAnalyzer* ptr = 0;
START_SECTION(ClusterAnalyzer())
{
	ptr = new ClusterAnalyzer();
	TEST_NOT_EQUAL(ptr, 0)
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

START_SECTION((std::vector< Real > averageSilhouetteWidth(std::vector< BinaryTreeNode > &tree, DistanceMatrix< Real > &original)))
{
	DistanceMatrix<Real> matrix(6,666);
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

	vector<Real> asw(5);
	asw[0]=0.170833f;
	asw[1]=0.309722f;
	asw[2]=0.306412f;
	asw[3]=0.125744f;
	asw[4]=0;

	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(2,3,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,1,0.6f));
	tree.push_back(BinaryTreeNode(0,1,0.7f));
	vector<Real> result = ptr->averageSilhouetteWidth(tree, matrix);
	TEST_EQUAL(result.size(), asw.size());
	for (Size i = 0; i < result.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.001);
			TEST_REAL_SIMILAR(result[i], asw[i]);
	}

}
END_SECTION

START_SECTION((std::vector< Real > dunnIndices(std::vector< BinaryTreeNode > &tree, DistanceMatrix< Real > &original, bool tree_from_singlelinkage=false)))
{
	DistanceMatrix<Real> matrix(6,666);
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
	tree.push_back(BinaryTreeNode(2,3,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,1,0.6f));
	tree.push_back(BinaryTreeNode(0,1,0.7f));
	vector<Real> di(5);
	di[0]=0.4f/0.3f; di[1]=0.5f/0.4f; di[2]=0.6f/0.8f; di[3]=0.7f/0.8f; di[4]=0.0f;
	vector<Real> result = ptr->dunnIndices(tree, matrix);
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

START_SECTION((void cut(Size cluster_quantity, std::vector< std::vector< Size > > &clusters, std::vector< BinaryTreeNode > &tree)))
{
	DistanceMatrix<Real> matrix(6,666);
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
	vector< vector<Size> > result;

	result.push_back(vector<Size>(a,a+3));
	result.push_back(vector<Size>(a+3,a+5));
	result.push_back(vector<Size>(a+5,a+6));

  vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(2,3,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,1,0.6f));
	tree.push_back(BinaryTreeNode(0,1,0.7f));
  ptr->cut(3, clusters, tree);
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

START_SECTION((std::vector< Real > cohesion(std::vector< std::vector< Size > > &clusters, DistanceMatrix< Real > &original)))
{
	DistanceMatrix<Real> matrix(6,666);
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

	vector<Real> cohesions;
	cohesions.push_back(0.533f);
	cohesions.push_back(0.4f);
	cohesions.push_back(0.7f);

	vector<Real> result = ptr->cohesion(clusters, matrix);
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

START_SECTION((Real averagePopulationAberration(Size cluster_quantity, std::vector< BinaryTreeNode > &tree)))
{
	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(2,3,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,1,0.6f));
	tree.push_back(BinaryTreeNode(0,1,0.7f));

	Real result = ptr->averagePopulationAberration(3, tree);
	TEST_REAL_SIMILAR(2.0/3.0, result);
}
END_SECTION

START_SECTION((String newickTree(std::vector< BinaryTreeNode > &tree, bool include_distance=false)))
{
	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(2,3,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,1,0.6f));
	tree.push_back(BinaryTreeNode(0,1,0.7f));

	String result = ptr->newickTree(tree);
	TEST_EQUAL(result,"( ( ( 0 , ( 1 , 2 ) ) , ( 3 , 4 ) ) , 5 )");
	result = ptr->newickTree(tree,true);
	TEST_EQUAL(result,"( ( ( 0:0.5 , ( 1:0.3 , 2:0.3 ):0.5 ):0.6 , ( 3:0.4 , 4:0.4 ):0.6 ):0.7 , 5:0.7 )");
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



