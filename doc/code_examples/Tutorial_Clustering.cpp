// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//

#include <OpenMS/ML/CLUSTERING/ClusterAnalyzer.h>
#include <OpenMS/ML/CLUSTERING/ClusterHierarchical.h>
#include <OpenMS/ML/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/ML/CLUSTERING/SingleLinkage.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <iostream>
#include <random>
#include <vector>

using namespace OpenMS;
using namespace std;



/// A functor, which provides a similarity value for two entities (here: doubles), in range [0, 1)
class LowLevelComparator
{
public:
  double operator()(const double first, const double second) const
  {
    // we just use a linear distance between them, i.e. the closer the values, the more similar they are
    auto distance = std::fabs(first - second);
    if (distance > 1) { throw Exception::InvalidRange(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION); }
    return 1 - distance;
  }
}; // end of LowLevelComparator

Int main()
{
  // data
  vector<double> data;
#if 1  // manual data
  data = {0.01, 0.02, 0.7, 0.3, 0.31};
#else  // random data
  const auto N = 5;
  std::mt19937 rng;                               // default constructed, seeded with fixed seed
  std::uniform_real_distribution<> dis(0.0, 1.0); // uniform values between [0, 1)
  std::generate_n(back_inserter(data), N, [&]() { return dis(rng); });
#endif

  // print raw data to console
  std::cout << "raw data: ";
  for_each(data.begin(), data.end(), [](auto elem) { std::cout << elem << ' '; });
  std::cout << '\n';
  // determines the distance between two data points
  LowLevelComparator llc;
  
  SingleLinkage sl; 
  // or try: 
  //CompleteLinkage sl;

  vector<BinaryTreeNode> tree;
  DistanceMatrix<float> dist; // will be filled
  ClusterHierarchical ch;
  ch.setThreshold(1); // maximal distance between clusters; default threshold = 1, i.e. full clustering
  // note: not all methods support a threshold, e.g. SingleLinkage requires t = 1.

  // do clustering.
  // Note: There are other overloads of this function for clustering spectra
  ch.cluster<double, LowLevelComparator>(data, llc, sl, tree, dist);

  // depending on the cluster method, the distance matrix may have shrunken, e.g. for complete linkage to the point where clustering was stopped
  std::cout << "distance matrix:\n" << dist << "\n\n"; 

  ClusterAnalyzer ca;
  std::cout << "binary tree in Newick format (numbers are indices into the data)";
  std::cout << ca.newickTree(tree) << std::endl;

  return 0;
} // end of main
