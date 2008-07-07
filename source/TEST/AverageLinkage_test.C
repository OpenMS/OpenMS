// -*- Mode: C++; tab-width: 2; -*-
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
#include <OpenMS/COMPARISON/CLUSTERING/AverageLinkage.h>
#include <OpenMS/DATASTRUCTURES/DistanceMatrix.h>
#include <vector>
//#include <iostream>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(AverageLinkage, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

AverageLinkage* ptr = 0;
CHECK(AverageLinkage())
{
	ptr = new AverageLinkage();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~AverageLinkage())
{
	delete ptr;
}
RESULT

CHECK((AverageLinkage(const AverageLinkage &source)))
{
	AverageLinkage al1;
	AverageLinkage copy(al1);
	TEST_EQUAL(copy.getName(), al1.getName());
}
RESULT

CHECK((AverageLinkage& operator=(const AverageLinkage &source)))
{
	AverageLinkage copy,al2;
	copy = al2;
	TEST_EQUAL(copy.getName(), al2.getName());
}
RESULT

CHECK((void cluster(const DistanceMatrix< double > &originalDist, DistanceMatrix< double > &actualDist, vector< vector< UInt > > &clusters, const String filepath="", const double threshold=1) const  throw (Exception::UnableToCreateFile, ClusterFunctor::InsufficientInput)))
{
	DistanceMatrix<double> matrix(5,666);
	matrix.setValue(0,1,1.0);
	matrix.setValue(0,2,2.8);
	matrix.setValue(0,3,3.6);
	matrix.setValue(0,4,4.2);
	matrix.setValue(1,2,2.2);
	matrix.setValue(1,3,2.8);
	matrix.setValue(1,4,3.6);
	matrix.setValue(2,3,1.0);
	matrix.setValue(2,4,1.4);
	matrix.setValue(3,4,1.0);
	DistanceMatrix<double> m(matrix);
	vector< vector<UInt> > result;
	UInt a[] = {0,1,2,3,4};
	vector< vector<UInt> > cluster;
	double th(3.1);
	AverageLinkage al3;
	al3.cluster(m,matrix,cluster,"",th);
	result.push_back(vector<UInt>(a,a+2));
	result.push_back(vector<UInt>(a+2,a+5));
	TEST_EQUAL(cluster.size(), result.size());
	for (UInt i = 0; i < cluster.size(); ++i)
	{
		TEST_EQUAL(cluster[i].size(), result[i].size());
		for (UInt j = 0; j < cluster[i].size(); ++j)
		{
			TEST_EQUAL(cluster[i][j], result[i][j]);
		}
	}
	th = 3.21;
	result.clear();
	result.push_back(vector<UInt>(a,a+5));
	DistanceMatrix<double> matrix2 = m;
	cluster.clear();
	AverageLinkage al4;
	al4.cluster(m,matrix2,cluster,"",th);
	TEST_EQUAL(cluster.size(), result.size());
	for (UInt i = 0; i < cluster.size(); ++i)
	{
		TEST_EQUAL(cluster[i].size(), result[i].size());
		for (UInt j = 0; j < cluster[i].size(); ++j)
		{
			TEST_EQUAL(cluster[i][j], result[i][j]);
		}
	}
	
}
RESULT

CHECK((static const String getProductName()))
{
	AverageLinkage al5;
	TEST_EQUAL(al5.getProductName(), "AverageLinkage")
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



