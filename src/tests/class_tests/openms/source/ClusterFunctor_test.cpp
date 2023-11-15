// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/AverageLinkage.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ClusterFunctor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(ClusterFunctor())
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(~ClusterFunctor())
{
  NOT_TESTABLE
}
END_SECTION


START_SECTION((ClusterFunctor(const ClusterFunctor &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((ClusterFunctor& operator=(const ClusterFunctor &source)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void operator()(DistanceMatrix< float > &original_distance, std::vector<BinaryTreeNode>& cluster_tree, const float threshold=1) const =0))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static void registerChildren()))
{
  ClusterFunctor* cfp = Factory<ClusterFunctor>::create("AverageLinkage");
  AverageLinkage* avl_nullPointer = nullptr;
  TEST_NOT_EQUAL( dynamic_cast<AverageLinkage*>(cfp) , avl_nullPointer)
  delete cfp;

  cfp = Factory<ClusterFunctor>::create("SingleLinkage");
  SingleLinkage* sl_nullPointer = nullptr;
  TEST_NOT_EQUAL( dynamic_cast<SingleLinkage*>(cfp) , sl_nullPointer)
  delete cfp;

  cfp = Factory<ClusterFunctor>::create("CompleteLinkage");
  CompleteLinkage* cl_nullPointer = nullptr;
  TEST_NOT_EQUAL( dynamic_cast<CompleteLinkage*>(cfp) , cl_nullPointer)
  delete cfp;
}
END_SECTION

START_SECTION(([ClusterFunctor::InsufficientInput] InsufficientInput(const char *file, int line, const char *function, const char *message="not enough data points to cluster anything")))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(([ClusterFunctor::InsufficientInput] virtual ~InsufficientInput()))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



