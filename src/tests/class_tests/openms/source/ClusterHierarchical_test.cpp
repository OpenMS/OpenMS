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
#include <OpenMS/ML/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/ML/CLUSTERING/SingleLinkage.h>
#include <OpenMS/KERNEL/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/BinnedSharedPeakCount.h>
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
 tree.push_back(BinaryTreeNode(0,1,0.00849858f));
 DistanceMatrix<float> matrix;

 ch.cluster(d,bspc,1.5,2,BinnedSpectrum::DEFAULT_BIN_OFFSET_LOWRES,sl,result, matrix);

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



