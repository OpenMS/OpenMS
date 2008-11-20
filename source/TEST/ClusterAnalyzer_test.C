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
// $Maintainer: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ClusterAnalyzer, "$Id$")

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

	vector<Real> asw(5);
	asw[0]=0.170833;
	asw[1]=0.309722;
	asw[2]=0.306412;
	asw[3]=0.125744;
	asw[4]=0;

	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3));
	tree.push_back(BinaryTreeNode(2,3,0.4));
	tree.push_back(BinaryTreeNode(0,1,0.5));
	tree.push_back(BinaryTreeNode(0,1,0.6));
	tree.push_back(BinaryTreeNode(0,1,0.7));
	vector<Real> result = ptr->averageSilhouetteWidth(tree, matrix);
	TEST_EQUAL(result.size(), asw.size());
	for (UInt i = 0; i < result.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.001);
			TEST_REAL_SIMILAR(result[i], asw[i]);
	}

}
END_SECTION

START_SECTION((std::vector< Real > dunnIndices(std::vector< BinaryTreeNode > &tree, DistanceMatrix< Real > &original, bool tree_from_singlelinkage)))
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

  vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3));
	tree.push_back(BinaryTreeNode(2,3,0.4));
	tree.push_back(BinaryTreeNode(0,1,0.5));
	tree.push_back(BinaryTreeNode(0,1,0.6));
	tree.push_back(BinaryTreeNode(0,1,0.7));
	vector<Real> di(5);
	di[0]=0.4/0.3; di[1]=0.5/0.4; di[2]=0.6/0.8; di[3]=0.7/0.8; di[4]=0.0;
	vector<Real> result = ptr->dunnIndices(tree, matrix);
	TEST_EQUAL(result.size(), di.size());
	for (UInt i = 0; i < result.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.001);
			TEST_REAL_SIMILAR(result[i], di[i]);
	}
	result = ptr->dunnIndices(tree, matrix, true);
	TEST_EQUAL(result.size(), di.size());
	for (UInt i = 0; i < result.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.001);
			TEST_REAL_SIMILAR(result[i], di[i]);
	}
}
END_SECTION

START_SECTION((void cut(std::vector< BinaryTreeNode > &tree, std::vector< std::vector< UInt > > &clusters, size_t cluster_quantity)))
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

  UInt a[] = {0,1,2,3,4,5};
	vector< vector<UInt> > clusters;
	vector< vector<UInt> > result;

	result.push_back(vector<UInt>(a,a+3));
	result.push_back(vector<UInt>(a+3,a+5));
	result.push_back(vector<UInt>(a+5,a+6));

  vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3));
	tree.push_back(BinaryTreeNode(2,3,0.4));
	tree.push_back(BinaryTreeNode(0,1,0.5));
	tree.push_back(BinaryTreeNode(0,1,0.6));
	tree.push_back(BinaryTreeNode(0,1,0.7));
  ptr->cut(3, clusters, tree);
	TEST_EQUAL(clusters.size(), result.size());
	for (UInt i = 0; i < clusters.size(); ++i)
	{
		TEST_EQUAL(clusters[i].size(), result[i].size());
		for (UInt j = 0; j < clusters[i].size(); ++j)
		{
			TEST_EQUAL(clusters[i][j], result[i][j]);
		}
	}
}
END_SECTION

START_SECTION((std::vector< Real > cohesion(std::vector< std::vector< UInt > > &clusters, DistanceMatrix< Real > &original)))
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

  UInt a[] = {0,1,2,3,4,5};
	vector< vector<UInt> > clusters;
	clusters.push_back(vector<UInt>(a,a+3));
	clusters.push_back(vector<UInt>(a+3,a+5));
	clusters.push_back(vector<UInt>(a+5,a+6));

	vector<Real> cohesions;
	cohesions.push_back((12.0+2.0/3.0)/100);
	cohesions.push_back(0);
	cohesions.push_back(0);

	vector<Real> result = ptr->cohesion(clusters, matrix);
	TEST_EQUAL(cohesions.size(), result.size());
	for (UInt i = 0; i < cohesions.size(); ++i)
	{
		TEST_EQUAL(cohesions[i], result[i]);
	}
}
END_SECTION

START_SECTION((Real averagePopulationAberration(size_t cluster_quantity, std::vector< BinaryTreeNode > &tree)))
{
	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3));
	tree.push_back(BinaryTreeNode(2,3,0.4));
	tree.push_back(BinaryTreeNode(0,1,0.5));
	tree.push_back(BinaryTreeNode(0,1,0.6));
	tree.push_back(BinaryTreeNode(0,1,0.7));

	Real result = ptr->averagePopulationAberration(3, tree);
	TEST_REAL_SIMILAR(2.0/3.0, result);
}
END_SECTION

START_SECTION((String newickTree(std::vector< BinaryTreeNode > &tree, bool include_distance)))
{
	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3));
	tree.push_back(BinaryTreeNode(2,3,0.4));
	tree.push_back(BinaryTreeNode(0,1,0.5));
	tree.push_back(BinaryTreeNode(0,1,0.6));
	tree.push_back(BinaryTreeNode(0,1,0.7));

	String result = ptr->newickTree(tree);
	TEST_EQUAL("( ( ( 0 , ( 1 , 2 ) ) , ( 3 , 4 ) ) , 5 )", result);
	result = ptr->newickTree(tree,true);
	TEST_EQUAL("( ( ( 0:0.5 , ( 1:0.3 , 2:0.3 ):0.5 ):0.6 , ( 3:0.4 , 4:0.4 ):0.6 ):0.7 , 5:0.7 )", result);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



