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
  namespace Interfaces { class IMSDataConsumer; }

  /// Dispatcher that detects whether an mzML contains spectra (SWATH/DIA)
  /// or chromatograms only (SRM/MRM) and forwards to the appropriate loader.
  struct OPENMS_DLLAPI TargetedDataFileLoader
  {
    /// Load a file and return Swath maps (may represent chromatogram-only SRM/MRM data)
    static std::vector<::OpenSwath::SwathMap> loadFile(const String& file,
                                                      const String& tmp,
                                                      std::shared_ptr<ExperimentalSettings>& exp_meta,
                                                      const String& readoptions,
                                                      Interfaces::IMSDataConsumer* plugin_consumer = nullptr);
  };

}
