// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
#include <OpenMS/CONCEPT/Factory.h>

using namespace std;

namespace OpenMS
{
  ClusterFunctor::ClusterFunctor() = default;

  ClusterFunctor::ClusterFunctor(const ClusterFunctor & /*source*/) = default;

  ClusterFunctor::~ClusterFunctor() = default;

  ClusterFunctor & ClusterFunctor::operator=(const ClusterFunctor & /*source*/) = default;

  void ClusterFunctor::registerChildren()
  {
    Factory<ClusterFunctor>::registerProduct(SingleLinkage::getProductName(), &SingleLinkage::create);
    Factory<ClusterFunctor>::registerProduct(CompleteLinkage::getProductName(), &CompleteLinkage::create);
    Factory<ClusterFunctor>::registerProduct(AverageLinkage::getProductName(), &AverageLinkage::create);
  }

  ClusterFunctor::InsufficientInput::InsufficientInput(const char * file, int line, const char * function, const char * message) throw() :
    BaseException(file, line, function, "ClusterFunctor::InsufficentInput", message)
  {
  }

  ClusterFunctor::InsufficientInput::~InsufficientInput() throw() = default;

}
