// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/ML/CLUSTERING/GridBasedClustering.h>

using namespace std;

namespace OpenMS
{
    
  MinimumDistance::MinimumDistance(const int &cluster_index, const int &nearest_neighbour_index, const double &distance)
  :cluster_index_(cluster_index), nearest_neighbour_index_(nearest_neighbour_index), distance_(distance)
  {
  }
  
  int MinimumDistance::getClusterIndex() const
  {
    return cluster_index_;
  }

  int MinimumDistance::getNearestNeighbourIndex() const
  {
    return  nearest_neighbour_index_;
  }

  bool MinimumDistance::operator<(const MinimumDistance& other) const
  {
    return distance_ < other.distance_;
  }
  
  bool MinimumDistance::operator>(const MinimumDistance& other) const
  {
    return distance_ > other.distance_;
  }
  
  bool MinimumDistance::operator==(const MinimumDistance& other) const
  {
    return distance_ == other.distance_;
  }

}
