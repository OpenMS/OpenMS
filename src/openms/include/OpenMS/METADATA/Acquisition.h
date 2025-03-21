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
      @brief Information about one raw data spectrum that was combined with several
      other raw data spectra.

      Although this class is basically a string value, it is needed to store important meta info for each raw data scan.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Acquisition :
    public MetaInfoInterface
  {
public:
    /// Constructor
    Acquisition() = default;
    /// Copy constructor
    Acquisition(const Acquisition &) = default;
    /// Move constructor
    Acquisition(Acquisition&&) = default;
    /// Destructor
    ~Acquisition() = default;

    /// Assignment operator
    Acquisition & operator=(const Acquisition &) = default;
    /// Move assignment operator
    Acquisition& operator=(Acquisition&&) & = default;

    /// Equality operator
    bool operator==(const Acquisition & rhs) const;
    /// Equality operator
    bool operator!=(const Acquisition & rhs) const;

    /// return the identifier/index/number of the acquisition
    const String & getIdentifier() const;
    /// sets the index/number of the scan
    void setIdentifier(const String & identifier);

protected:
    String identifier_;

  };
} // namespace OpenMS

