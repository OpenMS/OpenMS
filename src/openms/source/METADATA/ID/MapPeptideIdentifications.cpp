// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/METADATA/ID/MapPeptideIdentifications.h>

#include <algorithm>

namespace OpenMS
{
  const PeptideIdentificationList& MapPeptideIdentifications::at(const std::vector<Size>& path) const
  {
    static const PeptideIdentificationList none;
    const auto found = features.find(path);
    return found == features.end() ? none : found->second;
  }

  std::vector<const PeptideIdentification*> MapPeptideIdentifications::peptides() const
  {
    std::vector<const PeptideIdentification*> result;
    for (const auto& [path, items] : features)
    {
      for (const auto& item : items)
      {
        result.push_back(&item);
      }
    }
    for (const auto& item : unassigned)
    {
      result.push_back(&item);
    }
    return result;
  }

  void MapPeptideIdentifications::eraseFeature(const std::vector<Size>& path)
  {
    // the paths of a feature and its subordinates start with its path, and follow each other
    auto end = features.lower_bound(path);
    while (end != features.end() && end->first.size() >= path.size() && std::equal(path.begin(), path.end(), end->first.begin()))
    {
      ++end;
    }
    features.erase(features.lower_bound(path), end);
  }
} // namespace OpenMS
