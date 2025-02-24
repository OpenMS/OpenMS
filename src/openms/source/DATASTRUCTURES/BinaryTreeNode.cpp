// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/BinaryTreeNode.h>

namespace OpenMS
{

  BinaryTreeNode::BinaryTreeNode(const Size i, const Size j, const float x) :
    left_child(i), right_child(j), distance(x)
  {
  }

  BinaryTreeNode::BinaryTreeNode(const BinaryTreeNode& source)  = default;

  BinaryTreeNode::~BinaryTreeNode() = default;

  BinaryTreeNode& BinaryTreeNode::operator=(const BinaryTreeNode& source)
  {
    if (this != &source)
    {
      left_child = source.left_child;
      right_child = source.right_child;
      distance = source.distance;
    }
    return *this;
  }

}
