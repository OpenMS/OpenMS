// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/ClusterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/SingleLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/CompleteLinkage.h>
#include <OpenMS/COMPARISON/CLUSTERING/AverageLinkage.h>

using namespace std;

namespace OpenMS
{
  ClusterFunctor::ClusterFunctor() = default;

  ClusterFunctor::ClusterFunctor(const ClusterFunctor & /*source*/) = default;

  ClusterFunctor::~ClusterFunctor() = default;

  ClusterFunctor & ClusterFunctor::operator=(const ClusterFunctor & /*source*/) = default;

  ClusterFunctor::InsufficientInput::InsufficientInput(const char * file, int line, const char * function, const char * message) throw() :
    BaseException(file, line, function, "ClusterFunctor::InsufficentInput", message)
  {
  }

  ClusterFunctor::InsufficientInput::~InsufficientInput() throw() = default;

}
