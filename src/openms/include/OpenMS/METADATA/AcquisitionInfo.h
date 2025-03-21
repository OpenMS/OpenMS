// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/Acquisition.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Description of the combination of raw data to a single spectrum

    Specification for combining raw scans ( Acquisition ) into a single spectrum.
    A list of acquisitions from the original raw file can be specified.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI AcquisitionInfo :
    private std::vector<Acquisition>,
    public MetaInfoInterface
  {
private:
    typedef std::vector<Acquisition> ContainerType;

public:
    /// Constructor
    AcquisitionInfo() = default;
    /// Copy constructor
    AcquisitionInfo(const AcquisitionInfo&) = default;
    /// Move constructor
    AcquisitionInfo(AcquisitionInfo&&) = default;
    /// Destructor
    ~AcquisitionInfo() = default;

    /// Assignment operator
    AcquisitionInfo& operator=(const AcquisitionInfo&) = default;
    /// Move assignment operator
    AcquisitionInfo& operator=(AcquisitionInfo&&) & = default;

    /// Equality operator
    bool operator==(const AcquisitionInfo& rhs) const;
    /// Equality operator
    bool operator!=(const AcquisitionInfo& rhs) const;

    /// returns the method of combination
    const String& getMethodOfCombination() const;
    /// sets the method of combination
    void setMethodOfCombination(const String& method_of_combination);

    ///@name Export methods from private base std::vector<Acquisition>
    //@{
    using ContainerType::operator[];
    using ContainerType::begin;
    using ContainerType::end;
    using ContainerType::size;
    using ContainerType::push_back;
    using ContainerType::empty;
    using ContainerType::back;
    using ContainerType::insert;
    using ContainerType::resize;

    using ContainerType::iterator;
    using ContainerType::const_iterator;
    //@}

protected:
    String method_of_combination_;

  };
} // namespace OpenMS

