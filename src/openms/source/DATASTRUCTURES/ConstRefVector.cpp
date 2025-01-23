// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/DATASTRUCTURES/ConstRefVector.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  ConstRefVector<FeatureMap > default_constrefvector;
  ConstRefVector<FeatureMap >::Iterator default_constrefvector_iterator;
  ConstRefVector<FeatureMap >::ConstIterator default_constrefvector_constiterator;
}
