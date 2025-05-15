// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{

  /** @brief Elements of a binary tree used to represent a hierarchical clustering process

          strict indexing/topology is assumed, i.e. node no. x represents clusteringstep no. x
          left_child and right_child are each the lowest indices to elements of the merged clusters, distance is the distance of the two children
  */
  class OPENMS_DLLAPI BinaryTreeNode
  {
public:
    /// constructor
    BinaryTreeNode(const Size i, const Size j, const float x);

    /// destructor
    ~BinaryTreeNode();

    /// copy constructor
    BinaryTreeNode(const BinaryTreeNode& source);

    /// assignment operator
    BinaryTreeNode& operator=(const BinaryTreeNode& source);

    Size left_child;
    Size right_child;
    float distance;

private:
    BinaryTreeNode();
  };

}

