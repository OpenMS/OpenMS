// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
      @brief Scan window description

      @ingroup Metadata
  */
  struct OPENMS_DLLAPI ScanWindow :
    public MetaInfoInterface
  {
    /// Default constructor
    ScanWindow() = default;
    /// Copy constructor
    ScanWindow(const ScanWindow &) = default;
    /// Move constructor
    ScanWindow(ScanWindow&&) = default;
    /// Destructor
    ~ScanWindow() {}

    /// Equality operator
    bool operator==(const ScanWindow & source) const;
    /// Equality operator
    bool operator!=(const ScanWindow & source) const;

    /// Assignment operator
    ScanWindow & operator=(const ScanWindow &) = default;
    /// Move assignment operator
    ScanWindow& operator=(ScanWindow&&) & = default;

    /// Begin of the window
    double begin = 0.0;
    /// End of the window
    double end = 0.0;
  };

} // namespace OpenMS

