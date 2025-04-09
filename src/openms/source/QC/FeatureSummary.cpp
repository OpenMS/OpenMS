// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------


#include <OpenMS/QC/FeatureSummary.h>

using namespace std;

namespace OpenMS
{
  FeatureSummary::Result FeatureSummary::compute(const FeatureMap& feature_map)
  {
    FeatureSummary::Result result;
    float sum_rt_deviations = 0;
    UInt rt_count = 0;
    result.feature_count = feature_map.size();
    for (const auto& f : feature_map)
    {
      if (f.metaValueExists("rt_deviation"))
      {
        sum_rt_deviations += (float)f.getMetaValue("rt_deviation");
        rt_count += 1;
      }
    }

    // calculate mean rt shift (sec)
    if (rt_count != 0)
    {
      result.rt_shift_mean = sum_rt_deviations / rt_count;
    }
    else
    {
      result.rt_shift_mean = 0;
    }

    return result;
  }

  bool FeatureSummary::Result::operator==(const Result& rhs) const
  {
    return feature_count == rhs.feature_count && rt_shift_mean == rhs.rt_shift_mean;
  }

  /// Returns the name of the metric
  const String& FeatureSummary::getName() const
  {
    return name_;
  }

  /// Returns required file input i.e. MzML.
  /// This is encoded as a bit in a Status object.
  QCBase::Status FeatureSummary::requirements() const
  {
    return QCBase::Status(QCBase::Requires::PREFDRFEAT);
  }
} // namespace OpenMS
