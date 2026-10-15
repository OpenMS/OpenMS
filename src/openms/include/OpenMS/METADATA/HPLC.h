// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/METADATA/Gradient.h>

namespace OpenMS
{
  /**
    @brief Representation of a HPLC experiment

    It contains the description of instrument, the settings and the gradient.

    @ingroup Metadata
  */
  class OPENMS_DLLAPI HPLC
  {
public:
    /// Constructor
    HPLC();
    /// Copy constructor
    HPLC(const HPLC &) = default;
    /// Move constructor
    HPLC(HPLC&&) = default;
    /// Destructor
    ~HPLC();

    /// Assignment operator
    HPLC & operator=(const HPLC &) = default;
    /// Move assignment operator
    HPLC& operator=(HPLC&&) & = default;

    /// Equality operator
    bool operator==(const HPLC & source) const;
    /// Equality operator
    bool operator!=(const HPLC & source) const;

    /// returns a const reference to the instrument name
    const std::string & getInstrument() const;
    /// sets the instrument name
    void setInstrument(const std::string & instrument);

    /// returns a const reference to the column description
    const std::string & getColumn() const;
    /// sets the column description
    void setColumn(const std::string & column);

    /// returns the temperature (in degree C)
    Int getTemperature() const;
    /// sets the temperature (in degree C)
    void setTemperature(Int temperature);

    /// returns the pressure (in bar)
    UInt getPressure() const;
    /// sets the pressure (in bar)
    void setPressure(UInt pressure);

    /// returns the flux (in microliter/sec)
    UInt getFlux() const;
    /// sets the flux (in microliter/sec)
    void setFlux(UInt flux);

    /// returns the comments
    std::string getComment() const;
    /// sets the comments
    void setComment(std::string comment);

    /// returns a const reference to the used gradient
    const Gradient & getGradient() const;
    /// returns a mutable reference to the used gradient
    Gradient & getGradient();
    /// sets the used gradient
    void setGradient(const Gradient & gradient);

protected:
    std::string instrument_;
    std::string column_;
    Int temperature_;
    Int pressure_;
    Int flux_;
    std::string comment_;
    Gradient gradient_;
  };

} // namespace OpenMS

