// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/ML/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/ML/CLUSTERING/SingleLinkage.h>
#include <OpenMS/ML/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/ML/CLUSTERING/AverageLinkage.h>
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



