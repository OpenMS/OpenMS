// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------
//

#include <OpenMS/ML/CLUSTERING/GridBasedCluster.h>

namespace OpenMS
{

GridBasedCluster::GridBasedCluster(const Point &centre, const Rectangle &bounding_box, const std::vector<int> &point_indices, const int &property_A, const std::vector<int> &properties_B)
: centre_(centre), bounding_box_(bounding_box), point_indices_(point_indices), property_A_(property_A), properties_B_(properties_B)
{
}

GridBasedCluster::GridBasedCluster(const Point &centre, const Rectangle &bounding_box, const std::vector<int> &point_indices)
: centre_(centre), bounding_box_(bounding_box), point_indices_(point_indices), property_A_(-1), properties_B_(point_indices.size(),-1)
{
}

const GridBasedCluster::Point& GridBasedCluster::getCentre() const
{
  return centre_;
}

const GridBasedCluster::Rectangle& GridBasedCluster::getBoundingBox() const
{
  return bounding_box_;
}

const std::vector<int>& GridBasedCluster::getPoints() const
{
  return point_indices_;
}

int GridBasedCluster::getPropertyA() const
{
  return property_A_;
}

const std::vector<int>& GridBasedCluster::getPropertiesB() const
{
  return properties_B_;
}

bool GridBasedCluster::operator<(const GridBasedCluster& other) const
{
  return centre_.getY() < other.centre_.getY();
}

bool GridBasedCluster::operator>(const GridBasedCluster& other) const
{
  return centre_.getY() > other.centre_.getY();
}

bool GridBasedCluster::operator==(const GridBasedCluster& other) const
{
  return centre_.getY() == other.centre_.getY();
}

}
