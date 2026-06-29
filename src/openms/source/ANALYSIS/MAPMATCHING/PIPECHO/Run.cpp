// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "Run.h"
#include "Util.h"

namespace OpenMS::PipEcho
{

/******************************************************************************/
bool Run::is_donor_feature(const Feature& feature)
{
  return Util::feature_hit(feature).has_value();
}

} // namespace OpenMS::PipEcho
