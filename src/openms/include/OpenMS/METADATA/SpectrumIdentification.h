// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/IdentificationHit.h>

namespace OpenMS
{
  /**
    @brief Represents a object which can store the information of an analysisXML instance

        //@todo docu (Andreas)

        @ingroup Metadata
  */
  class OPENMS_DLLAPI SpectrumIdentification :
    public MetaInfoInterface
  {
public:

    /// @name constructors,destructors,assignment operator
    //@{
    /// Default constructor
    SpectrumIdentification() = default;
    /// Destructor
    virtual ~SpectrumIdentification();
    /// Copy constructor
    SpectrumIdentification(const SpectrumIdentification &) = default;
    /// Move constructor
    SpectrumIdentification(SpectrumIdentification&&) = default;
    /// Assignment operator
    SpectrumIdentification & operator=(const SpectrumIdentification &) = default;
    /// Move assignment operator
    SpectrumIdentification& operator=(SpectrumIdentification&&) & = default;
    /// Equality operator
    bool operator==(const SpectrumIdentification & rhs) const;
    /// Inequality operator
    bool operator!=(const SpectrumIdentification & rhs) const;
    //@}

    // @name Accessors
    //@{
    /// sets the identification hits of this spectrum identification (corresponds to single peptide hit in the list)
    void setHits(const std::vector<IdentificationHit> & hits);

    /// adds a single identification hit to the hits
    void addHit(const IdentificationHit & hit);

    /// returns the identification hits of this spectrum identification
    const std::vector<IdentificationHit> & getHits() const;
    //@}

protected:

    String id_; ///< Identifier
    std::vector<IdentificationHit> hits_; ///< Single peptide hits
  };

} //namespace OpenMS
