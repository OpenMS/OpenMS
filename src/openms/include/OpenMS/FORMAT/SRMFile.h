// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <memory>

namespace OpenMS
{
  /// Minimal SRM file loader returning a single SwathMap wrapping the chromatogram container
  struct SRMFile
  {
    static std::vector<::OpenSwath::SwathMap> loadMzML(const String& file,
                                                    const String& tmp,
                                                    std::shared_ptr<ExperimentalSettings>& exp_meta);
  };

}
