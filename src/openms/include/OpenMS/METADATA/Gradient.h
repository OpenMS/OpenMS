// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Representation of a HPLC gradient

    It consists of several eluents and timepoints.
    Linear behaviour between timepoints is assumed.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI Gradient
  {
public:
    /// Constructor
    Gradient() = default;
    /// Copy constructor
    Gradient(const Gradient &) = default;
    /// Move constructor
    Gradient(Gradient&&) = default;
    /// Destructor
    ~Gradient();

    /// Assignment operator
    Gradient & operator=(const Gradient &) = default;
    /// Move assignment operator
    Gradient& operator=(Gradient&&) & = default;

    /// Equality operator
    bool operator==(const Gradient & source) const;
    /// Equality operator
    bool operator!=(const Gradient & source) const;

    /**
        @brief Adds an eluent at the end of the eluent array

        @exception Exception::InvalidValue is thrown if the same eluent name is used twice.
    */
    void addEluent(const String & eluent);
    /// removes all eluents
    void clearEluents();
    /// returns a const reference to the list of eluents
    const std::vector<String> & getEluents() const;

    /**
        @brief Adds a timepoint at the end of the timepoint array

        @exception Exception::OutOfRange is thrown if the new timepoint is before the last timepoint.
    */
    void addTimepoint(Int timepoint);
    /// removes all timepoints
    void clearTimepoints();
    /// returns a const reference to the list of timepoints
    const std::vector<Int> & getTimepoints() const;

    /**
        @brief sets the percentage of eluent @p eluent at timepoint @p timepoint

        @exception Exception::InvalidValue is thrown if the eluent, timepoint or percentage is invalid.
    */
    void setPercentage(const String & eluent, Int timepoint, UInt percentage);

    /**
        @brief returns a const reference to the percentages

        First dimension of the vector is the eluents, second dimension is the timepoints.
    */
    const std::vector<std::vector<UInt> > & getPercentages() const;

    /**
        @brief returns the percentage of an @p eluent at a @p timepoint

        @exception Exception::InvalidValue is thrown if the eluent or timepoint is invalid.
    */
    UInt getPercentage(const String & eluent, Int timepoint) const;

    /// sets all percentage values to 0
    void clearPercentages();

    /// checks if the percentages of all timepoints add up to 100%
    bool isValid() const;

protected:
    std::vector<String> eluents_;
    std::vector<Int> times_;
    // first dimension is eluents, second is times
    std::vector<std::vector<UInt> > percentages_;
  };

} // namespace OpenMS

