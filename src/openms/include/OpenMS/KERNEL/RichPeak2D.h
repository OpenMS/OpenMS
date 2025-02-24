// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/CONCEPT/UniqueIdInterface.h>

namespace OpenMS
{

  /**
    @brief A 2-dimensional raw data point or peak with meta information.

    This data structure is intended for continuous data or peak data.
    If you do not need to annotated single peaks with meta data, use Peak2D instead.

    @ingroup Kernel
  */
  class OPENMS_DLLAPI RichPeak2D :
    public Peak2D,
    public MetaInfoInterface,
    public UniqueIdInterface
  {
public:

    /// Default constructor
    RichPeak2D() :
      Peak2D(),
      MetaInfoInterface(),
      UniqueIdInterface()
    {}

    /// Copy constructor
    RichPeak2D(const RichPeak2D& p) = default;

    /// Constructor from Peak2D
    explicit RichPeak2D(const Peak2D& p) :
      Peak2D(p),
      MetaInfoInterface()
    {
      UniqueIdInterface::clearUniqueId();
    }

    /// Member constructor
    explicit RichPeak2D(const PositionType& pos, const IntensityType in) :
      Peak2D(pos, in),
      MetaInfoInterface()
    {}

    /// Move constructor
    RichPeak2D(RichPeak2D&& p) = default;

    /// Destructor
    ~RichPeak2D() override
    {}

    /// Assignment operator
    RichPeak2D & operator=(const RichPeak2D& rhs) = default;

    /// Move Assignment operator
    RichPeak2D & operator=(RichPeak2D&& rhs) & = default;

    /// Assignment operator
    RichPeak2D & operator=(const Peak2D& rhs)
    {
      if (this == &rhs) return *this;

      Peak2D::operator=(rhs);
      clearMetaInfo();
      UniqueIdInterface::clearUniqueId();

      return *this;
    }

    /// Equality operator
    bool operator==(const RichPeak2D& rhs) const
    {
      return Peak2D::operator==(rhs) &&
             MetaInfoInterface::operator==(rhs) &&
             UniqueIdInterface::operator==(rhs);
    }

    /// Equality operator
    bool operator!=(const RichPeak2D& rhs) const
    {
      return !(operator==(rhs));
    }

  };

} // namespace OpenMS

