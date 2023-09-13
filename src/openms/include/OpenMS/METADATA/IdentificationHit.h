// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>

namespace OpenMS
{
  /**
    @brief Represents a object which can store the information of an analysisXML instance

    //@todo docu (Andreas)

    @ingroup Metadata
  */
  class OPENMS_DLLAPI IdentificationHit :
    public MetaInfoInterface
  {
public:

    /// @name constructors,destructors,assignment operator
    //@{

    /// Default constructor
    IdentificationHit();
    /// Copy constructor
    IdentificationHit(const IdentificationHit &) = default;
    /// Destructor
    virtual ~IdentificationHit();
    /// Move constructor
    IdentificationHit(IdentificationHit&&) = default;

    /// Assignment operator
    IdentificationHit & operator=(const IdentificationHit &) = default;
    /// Move assignment operator
    IdentificationHit& operator=(IdentificationHit&&) & = default;

    /// Equality operator
    bool operator==(const IdentificationHit & rhs) const;
    /// Inequality operator
    bool operator!=(const IdentificationHit & rhs) const;
    //@}

    /// @name Accessors
    //@{
    /// sets the identifier
    void setId(const String & id);

    /// returns the id
    const String & getId() const;

    /// sets the charge state of the peptide
    void setCharge(Int charge);

    /// returns the charge state
    Int getCharge() const;

    /// sets the calculated mass to charge ratio
    void setCalculatedMassToCharge(double mz);

    /// returns the calculated mass to charge ratio
    double getCalculatedMassToCharge() const;

    /// sets the experimental mass to charge ratio
    void setExperimentalMassToCharge(double mz);

    /// returns the experimental mass to charge
    double getExperimentalMassToCharge() const;

    /// sets the name
    void setName(const String & name);

    /// returns the name
    const String & getName() const;

    /// sets whether the peptide passed the threshold
    void setPassThreshold(bool pass);

    /// returns whether the peptide passed the threshold
    bool getPassThreshold() const;

    /// set the rank of the peptide
    void setRank(Int rank);

    /// returns the rank of the peptide
    Int getRank() const;
    //@}


protected:

    String id_;                                     ///< identifier
    Int charge_;                                    ///< peptide charge
    double calculated_mass_to_charge_;         ///< calculated mass to charge ratio
    double experimental_mass_to_charge_;         ///< experimental mass to charge ratio
    String name_;                               ///< name
    bool pass_threshold_;               ///< pass threshold
    Int rank_;                                      ///< rank of the peptide
  };

} //namespace OpenMS

