// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/ConversionHelper.h>

namespace OpenMS
{
  void MapConversion::convert(UInt64 const input_map_index,
                              PeakMap& input_map,
                              ConsensusMap& output_map,
                              Size n)
  {
    output_map.clear(true);

    // see @todo above
    output_map.setUniqueId();

    input_map.updateRanges(1);
    if (n > input_map.getSize())
    {
      n = input_map.getSize();
    }
    output_map.reserve(n);
    std::vector<Peak2D> tmp;
    tmp.reserve(input_map.getSize());

    // TODO Avoid tripling the memory consumption by this call
    input_map.get2DData(tmp);

    std::partial_sort(tmp.begin(),
                      tmp.begin() + n,
                      tmp.end(),
                      [](auto &left, auto &right) {Peak2D::IntensityLess cmp; return cmp(right, left);});

    for (Size element_index = 0; element_index < n; ++element_index)
    {
      output_map.push_back(ConsensusFeature(input_map_index,
                                            tmp[element_index],
                                            element_index));
    }

    output_map.getColumnHeaders()[input_map_index].size = n;
    output_map.updateRanges();
  }

  void MapConversion::convert(ConsensusMap const& input_map,
                              const bool keep_uids,
                              FeatureMap& output_map)
  {
    output_map.clear(true);
    output_map.resize(input_map.size());
    output_map.DocumentIdentifier::operator=(input_map);

    if (keep_uids)
    {
      output_map.UniqueIdInterface::operator=(input_map);
    }
    else
    {
      output_map.setUniqueId();
    }
    output_map.setProteinIdentifications(input_map.getProteinIdentifications());
    output_map.setUnassignedPeptideIdentifications(input_map.getUnassignedPeptideIdentifications());

    for (Size i = 0; i < input_map.size(); ++i)
    {
      Feature& f = output_map[i];
      const ConsensusFeature& c = input_map[i];
      f.BaseFeature::operator=(c);
      if (!keep_uids)
      {
        f.setUniqueId();
      }
    }

    output_map.updateRanges();
  }

  void MapConversion::convert(UInt64 const input_map_index,
                              FeatureMap const& input_map,
                              ConsensusMap& output_map,
                              Size n)
  {
    if (n > input_map.size())
    {
      n = input_map.size();
    }

    output_map.clear(true);
    output_map.reserve(n);

    // An arguable design decision, see above.
    output_map.setUniqueId(input_map.getUniqueId());

    for (UInt64 element_index = 0; element_index < n; ++element_index)
    {
      output_map.push_back(ConsensusFeature(input_map_index, input_map[element_index]));
    }
    output_map.getColumnHeaders()[input_map_index].size = static_cast<Size>(input_map.size());
    output_map.setProteinIdentifications(input_map.getProteinIdentifications());
    output_map.setUnassignedPeptideIdentifications(input_map.getUnassignedPeptideIdentifications());
    output_map.updateRanges();
  }

} // namespace OpenMS
