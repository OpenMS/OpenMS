// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MS1LabelState.h>

namespace OpenMS::MS1LabelState
{
  const std::vector<std::string>& keys()
  {
    static const std::vector<std::string> k{LABELED_SEQUENCE, REMOVED_LABELS, LABEL_CHANNEL};
    return k;
  }

  bool hasMatchedSequence(const PeptideHit& hit)
  {
    return hit.metaValueExists(LABELED_SEQUENCE);
  }

  AASequence matchedSequence(const PeptideHit& hit)
  {
    if (!hit.metaValueExists(LABELED_SEQUENCE)) { return hit.getSequence(); }
    return AASequence::fromString(hit.getMetaValue(LABELED_SEQUENCE).toString());
  }

  PeptideHit withMatchedSequence(const PeptideHit& hit)
  {
    PeptideHit matched = hit;
    if (hit.metaValueExists(LABELED_SEQUENCE))
    {
      matched.setSequence(AASequence::fromString(hit.getMetaValue(LABELED_SEQUENCE).toString()));
    }
    return matched;
  }
} // namespace OpenMS::MS1LabelState
