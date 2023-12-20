// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/METADATA/SpectrumIdentification.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Represents a object which can store the information of an analysisXML instance

        //@todo docu (Andreas)

        @ingroup Metadata
  */
  class OPENMS_DLLAPI Identification :
    public MetaInfoInterface
  {
public:

    /// @name constructors,destructors,assignment operator
    //@{

    /// Default constructor
    Identification() = default;
    /// Copy constructor
    Identification(const Identification & source) = default;
    /// Move constructor
    Identification(Identification&&) = default;
    /// Destructor
    virtual ~Identification();

    /// Assignment operator
    Identification & operator=(const Identification & source) = default;
    /// Move assignment operator
    Identification& operator=(Identification&&) & = default;

    /// Equality operator
    bool operator==(const Identification & rhs) const;
    /// Inequality operator
    bool operator!=(const Identification & rhs) const;
    //@}

    /// @name Accessors
    //@{
    /// sets the date and time the file was written
    void setCreationDate(const DateTime & date);

    /// returns the date and time the file was created
    const DateTime & getCreationDate() const;

    /// sets the spectrum identifications
    void setSpectrumIdentifications(const std::vector<SpectrumIdentification> & ids);

    /// adds a spectrum identification
    void addSpectrumIdentification(const SpectrumIdentification & id);

    /// returns the spectrum identifications stored
    const std::vector<SpectrumIdentification> & getSpectrumIdentifications() const;
    //@}
protected:

    String id_; ///< Identifier
    DateTime creation_date_; ///< Date and time the search was performed
    std::vector<SpectrumIdentification> spectrum_identifications_;

  };

} //namespace OpenMS

