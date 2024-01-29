// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/COMPARISON/CLUSTERING/ClusterHierarchical.h>

//using namespace std;

namespace OpenMS
{
  UnnormalizedComparator::UnnormalizedComparator(const char * file, int line, const char * function, const char * message) throw() :
    BaseException(file, line, function, "ClusterHierarchical::UnnormalizedComparator", message)
  {
  }

  UnnormalizedComparator::~UnnormalizedComparator() throw() = default;

}
