// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/COMPARISON/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace OpenMS;
using namespace std;


class LowLevelComparator
{
public:
  double operator()(const double first, const double second) const
  {
    double x, y;
    x = min(second, first);
    y = max(first, second);
    if ((y - x) > 1)
    {
      throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }
    return 1 - (y - x);
  }

}; // end of LowLevelComparator

Int main()
{
  // data
  vector<double> data; // must be filled

  srand(333);
  Size nr = 12;
  data.resize(nr);
  for (Size i = 0; i < nr; ++i)
  {
    data[i] = (double)rand() / RAND_MAX;
  }

  LowLevelComparator llc;
  CompleteLinkage sl;
  vector<BinaryTreeNode> tree;
  DistanceMatrix<float> dist; // will be filled
  ClusterHierarchical ch;
  ch.setThreshold(0.15);

  // clustering
  ch.cluster<double, LowLevelComparator>(data, llc, sl, tree, dist);

  ClusterAnalyzer ca;
  std::cout << ca.newickTree(tree) << std::endl;

  return 0;
} //end of main
