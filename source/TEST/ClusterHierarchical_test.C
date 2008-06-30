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
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <vector>
#include <algorithm>
///////////////////////////

using namespace OpenMS;
using namespace std;


class lowlevelComparator
{
	public:
	double operator()(const UInt first, const UInt second) const
	{
		UInt x,y;
		x = min(second,first);
		y = max(first,second);
		
		switch (x)
		{
			case 0: 
				switch(y)
				{
					default:
						return 666; 
						break;
					case 1:
						return 0; 
						break;
					case 2:
						return -1.8; 
						break;
					case 3:
						return -2.6; 
						break;
					case 4:
						return -3.2; 
						break;
				} 
			break;
			case 1:
				switch(y)
				{
					default:
						return 666; 
						break;
					case 2:
						return -1.2; 
						break;
					case 3:
						return -1.8; 
						break;
					case 4:
						return -2.6; 
						break;
				} 
			
			break;
			case 2:
				switch(y)
				{
					default:
						return 666; 
						break;
					case 3:
						return 0; 
						break;
					case 4:
						return -0.4; 
						break;
				} 
			
			break;
			case 3:
				switch(y)
				{
					default:
						return 666; 
						break;
					case 4:
						return 0; 
						break;
				} 
			
			break;
			default: 
				return 666;
				break;
		}
	}	
};


START_TEST(ClusterHierarchical, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ClusterHierarchical* ptr = 0;
CHECK(ClusterHierarchical())
{
	ptr = new ClusterHierarchical();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~ClusterHierarchical())
{
	delete ptr;
}
RESULT

CHECK((ClusterHierarchical(double x=1.0)))
{
  ClusterHierarchical ch(66.6);
	TEST_EQUAL(ch.getThreshold(),66.6)
}
RESULT

CHECK((ClusterHierarchical(const ClusterHierarchical &source)))
{
	ClusterHierarchical copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getThreshold(), ptr->getThreshold());
}
RESULT

CHECK((double getThreshold()))
{
	ptr->setThreshold(0.666);
	TEST_EQUAL(ptr->getThreshold(),0.666);
}
RESULT

CHECK((void setThreshold(double x)))
{
	ptr->setThreshold(0.666);
	TEST_EQUAL(ptr->getThreshold(),0.666);
}
RESULT
	
UInt a[] = {0,1,2,3,4};
vector<UInt> d(a,a+5);
lowlevelComparator lc;
SingleLinkage sl;
vector< vector<UInt> > result;
result.push_back(vector<UInt>(a,a+2));
result.push_back(vector<UInt>(a+2,a+5));

CHECK((template <typename Data, typename SimilarityComparator> void clusterForVector(vector< Data > &data, const SimilarityComparator &comparator, const ClusterFunctor &clusterer, vector< vector< UInt > > &clusters)))
{
	vector< vector<UInt> > r;
	ptr->setThreshold(2.2);
	ptr->clusterForVector<UInt,lowlevelComparator>(d,lc,sl,r);
	TEST_EQUAL(r.size(), result.size());
	for (UInt i = 0; i < r.size(); ++i)
	{
			TEST_EQUAL(r[i].size(), result[i].size());
			for (UInt j = 0; j < r[i].size(); ++j)
			{
				TEST_EQUAL(r[i][j], result[i][j]);
			}
	}
}
RESULT

CHECK((template <typename Data, typename SimilarityComparator> void clusterForDendrogramm(const vector< Data > &data, const SimilarityComparator &comparator, const ClusterFunctor &clusterer, vector< vector< UInt > > &clusters, const String &filepath)))
{
	vector< vector<UInt> > r;
	String s;
	NEW_TMP_FILE(s);
	const String string(s);
	ptr->setThreshold(2.2);
	ptr->clusterForDendrogramm<UInt,lowlevelComparator>(d,lc,sl,r,string);
	TEST_EQUAL(r.size(), result.size());
	for (UInt i = 0; i < r.size(); ++i)
	{
			TEST_EQUAL(r[i].size(), result[i].size());
			for (UInt j = 0; j < r[i].size(); ++j)
			{
				TEST_EQUAL(r[i][j], result[i][j]);
			}
	}
}
RESULT

CHECK((static const String getName()))
{
	TEST_EQUAL(ptr->getName(), "ClusterHierarchical");
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



