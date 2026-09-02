// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/MS1LabelState.h>

#include <OpenMS/METADATA/MetaInfo.h>
#include <OpenMS/METADATA/MetaInfoRegistry.h>

namespace OpenMS::MS1LabelState
{
  const std::vector<std::string>& keys()
  {
    static const std::vector<std::string> k{LABELED_SEQUENCE, REMOVED_LABELS, LABEL_CHANNEL};
    return k;
  }

  Keys::Keys() :
    labeled_sequence(MetaInfo::registry().getIndex(LABELED_SEQUENCE)),
    removed_labels(MetaInfo::registry().getIndex(REMOVED_LABELS)),
    label_channel(MetaInfo::registry().getIndex(LABEL_CHANNEL))
  {
  }

  bool hasMatchedSequence(const PeptideHit& hit)
  {
    return hit.metaValueExists(LABELED_SEQUENCE);
  }

  bool hasMatchedSequence(const PeptideHit& hit, const Keys& keys)
  {
    // UInt(-1) is guarded explicitly: the index overload of metaValueExists() does not special-case
    // the "not registered" sentinel the way the string overload does.
    return keys.labeled_sequence != static_cast<UInt>(-1) && hit.metaValueExists(keys.labeled_sequence);
  }

  AASequence matchedSequence(const PeptideHit& hit, const Keys& keys)
  {
    if (!hasMatchedSequence(hit, keys)) { return hit.getSequence(); }
    return AASequence::fromString(hit.getMetaValue(keys.labeled_sequence).toString());
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
