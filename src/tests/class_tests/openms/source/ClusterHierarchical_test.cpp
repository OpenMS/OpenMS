// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSharedPeakCount.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/DTAFile.h>

#include <vector>
#include <algorithm>
///////////////////////////

using namespace OpenMS;
using namespace std;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunreachable-code"

class LowlevelComparator
{
	public:
	double operator()(const Size first, const Size second) const
	{
		Size x,y;
		x = min(second,first);
		y = max(first,second);

		switch (x)
		{
			case 0:
				switch(y)
				{
					default:
						return 0;
						break;
					case 1:
						return 1-0.5;
						break;
					case 2:
						return 1-0.8;
						break;
					case 3:
						return 1-0.6;
						break;
					case 4:
						return 1-0.8;
						break;
					case 5:
						return 1-0.7;
						break;
				}
			break;
			case 1:
				switch(y)
				{
					default:
						return 0;
						break;
					case 2:
						return 1-0.3;
						break;
					case 3:
						return 1-0.8;
						break;
					case 4:
						return 1-0.8;
						break;
					case 5:
						return 1-0.8;
						break;
				}

			break;
			case 2:
				switch(y)
				{
					default:
						return 0;
						break;
					case 3:
						return 1-0.8;
						break;
					case 4:
						return 1-0.8;
						break;
					case 5:
						return 1-0.8;
						break;
				}

			break;
			case 3:
				switch(y)
				{
					default:
						return 0;
						break;
					case 4:
						return 1-0.4;
						break;
					case 5:
						return 1-0.8;
						break;
				}

			break;
			case 4:
				switch(y)
				{
					default:
						return 0;
						break;
					case 5:
						return 1-0.8;
						break;
				}

			break;
			default:
				return 666;
				break;

		}
	}
};

#pragma clang diagnostic pop

START_TEST(ClusterHierarchical, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ClusterHierarchical* ptr = nullptr;
ClusterHierarchical* nullPointer = nullptr;
START_SECTION(ClusterHierarchical())
{
	ptr = new ClusterHierarchical();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ClusterHierarchical())
{
	delete ptr;
}
END_SECTION

START_SECTION((ClusterHierarchical(const ClusterHierarchical &source)))
{
	ClusterHierarchical ch;
	ch.setThreshold(66.6);
	ClusterHierarchical copy(ch);
	TEST_EQUAL(copy.getThreshold(), 66.6);
}
END_SECTION

START_SECTION((double getThreshold()))
{
	ClusterHierarchical ch;
	ch.setThreshold(0.666);
	TEST_EQUAL(ch.getThreshold(),0.666);
}
END_SECTION

START_SECTION((void setThreshold(double x)))
{
	ClusterHierarchical ch;
	ch.setThreshold(0.666);
	TEST_EQUAL(ch.getThreshold(),0.666);
}
END_SECTION

START_SECTION((template <typename Data, typename SimilarityComparator> void cluster(std::vector< Data > &data, const SimilarityComparator &comparator, const ClusterFunctor &clusterer, std::vector<BinaryTreeNode>& cluster_tree, DistanceMatrix<float>& original_distance)))
{
	vector<Size> d(6,0);
	for (Size i = 0; i<d.size(); ++i)
	{
		d[i]=i;
	}
	ClusterHierarchical ch;
	LowlevelComparator lc;
	SingleLinkage sl;
	vector< BinaryTreeNode > result;
	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.3f));
	tree.push_back(BinaryTreeNode(3,4,0.4f));
	tree.push_back(BinaryTreeNode(0,1,0.5f));
	tree.push_back(BinaryTreeNode(0,3,0.6f));
	tree.push_back(BinaryTreeNode(0,5,0.7f));
	DistanceMatrix<float> matrix;

	ch.cluster<Size,LowlevelComparator>(d,lc,sl,result, matrix);

	TEST_EQUAL(tree.size(), result.size());
	for (Size i = 0; i < tree.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_EQUAL(tree[i].left_child, result[i].left_child);
			TEST_EQUAL(tree[i].right_child, result[i].right_child);
			TEST_REAL_SIMILAR(tree[i].distance, result[i].distance);
	}
}
END_SECTION

START_SECTION((void cluster(std::vector<PeakSpectrum>& data, const BinnedSpectrumCompareFunctor& comparator, double sz, UInt sp, const ClusterFunctor& clusterer, std::vector<BinaryTreeNode>& cluster_tree, DistanceMatrix<float>& original_distance)))
{

	PeakSpectrum s1, s2, s3;
	Peak1D peak;

	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
	s2 = s1;
	s3 = s1;
	s2.pop_back();
	s3.pop_back();
	peak.setMZ(666.66);
	peak.setIntensity(999.99f);
	s2.push_back(peak);
	s2.sortByPosition();
	s3.push_back(peak);
	s3.sortByPosition();

	vector<PeakSpectrum> d(3);
	d[0] = s1; d[1] = s2; d[2] = s3;
	ClusterHierarchical ch;
	BinnedSharedPeakCount bspc;
	SingleLinkage sl;
	vector< BinaryTreeNode > result;
	vector< BinaryTreeNode > tree;
	tree.push_back(BinaryTreeNode(1,2,0.0));
	tree.push_back(BinaryTreeNode(0,1,0.0086f));
	DistanceMatrix<float> matrix;

	ch.cluster(d,bspc,1.5,2,sl,result, matrix);

	TEST_EQUAL(tree.size(), result.size());
	for (Size i = 0; i < tree.size(); ++i)
	{
			TOLERANCE_ABSOLUTE(0.0001);
			TEST_EQUAL(tree[i].left_child, result[i].left_child);
			TEST_EQUAL(tree[i].right_child, result[i].right_child);
			TEST_REAL_SIMILAR(tree[i].distance, result[i].distance);
	}
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



