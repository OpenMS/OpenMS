// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------


#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/QC/SpectrumCount.h>
using namespace std;

namespace OpenMS
{

  map<Size, UInt> SpectrumCount::compute(const MSExperiment& exp)
  {
    map<Size, UInt> counts;
    for (const auto& spectrum : exp)
    {
      const Size level = spectrum.getMSLevel();
      ++counts[level]; // count MS level
    }
    return counts;
  }

  /// Returns the name of the metric
  const String& SpectrumCount::getName() const
  {
    return name_;
  }

  /// Returns required file input i.e. MzML.
  /// This is encoded as a bit in a Status object.
  QCBase::Status SpectrumCount::requirements() const
  {
    return QCBase::Status(QCBase::Requires::RAWMZML);
  }
} // namespace OpenMS
