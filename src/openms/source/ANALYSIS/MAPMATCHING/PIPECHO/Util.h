// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include <optional>

#include "OpenMS/KERNEL/Feature.h"
#include "OpenMS/METADATA/PeptideHit.h"

namespace OpenMS {
namespace PipEcho {
namespace Util {

  /****************************************************************************/
  /**
   * Return the first peptide hit from a feature.
   *
   * This code is here to isolate checking hits because the process
   * might change in a future version of OpenMS.
   */
 std::optional<PeptideHit> feature_hit(const Feature&);

}}} // Name spaces
