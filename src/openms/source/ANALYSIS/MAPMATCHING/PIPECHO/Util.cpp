// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "Util.h"

namespace OpenMS {
namespace PipEcho {
namespace Util {

  /****************************************************************************/
  std::optional<PeptideHit> feature_hit(const Feature& feature) {
    const PeptideIdentificationList& peps = feature.getPeptideIdentifications();

    // FIXME: Is this the correct way to know if a feature has an
    // ID?  There seems to be two ways to store an ID on a feature,
    // but most code I see uses the "old way".
    if (!peps.empty() && !peps[0].getHits().empty()) {
      return peps[0].getHits()[0];
    }

    return {};
  }

}}} // Name spaces
