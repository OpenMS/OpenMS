// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/CVTermList.h>

namespace OpenMS
{
  /**
      @brief Product meta information.

      This class describes the product isolation window for special scan types, such as MRM.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Product :
    public CVTermList
  {

public:

    /// Constructor
    Product() = default;
    /// Copy constructor
    Product(const Product &) = default;
    /// Move constructor
    Product(Product&&) = default;
    /// Destructor
    ~Product() override = default;

    /// Assignment operator
    Product & operator=(const Product &) = default;
    /// Move assignment operator
    Product& operator=(Product&&) & = default;

    /// Equality operator
    bool operator==(const Product & rhs) const;
    /// Equality operator
    bool operator!=(const Product & rhs) const;

    /// returns the target m/z
    double getMZ() const;
    /// sets the target m/z
    void setMZ(double mz);

    /// returns the lower offset from the target m/z
    double getIsolationWindowLowerOffset() const;
    /// sets the lower offset from the target m/z
    void setIsolationWindowLowerOffset(double bound);

    /// returns the upper offset from the target m/z
    double getIsolationWindowUpperOffset() const;
    /// sets the upper offset from the target m/z
    void setIsolationWindowUpperOffset(double bound);

protected:

    double mz_ = 0.0;
    double window_low_ = 0.0;
    double window_up_ = 0.0;
  };
} // namespace OpenMS

