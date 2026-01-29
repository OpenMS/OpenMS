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
  /**
    @brief Minimal SRM/MRM file loader returning a single SwathMap wrapping the chromatogram container.
  */
  struct OPENMS_DLLAPI MRMFile
  {
    /**
      @brief Load an SRM/MRM mzML file and return a single SwathMap.

      @param[in] file Input mzML file.
      @param[in] tmp Temporary directory (for cached data).
      @param[out] exp_meta Experimental settings extracted from the file.
      @return Swath maps (single-entry for SRM).
    */
    static std::vector<::OpenSwath::SwathMap> loadMzML(const String& file,
                                                    const String& tmp,
                                                    std::shared_ptr<ExperimentalSettings>& exp_meta);
  };

}
