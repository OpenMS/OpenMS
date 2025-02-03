// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/SampleTreatment.h>

#include <OpenMS/DATASTRUCTURES/String.h>

namespace OpenMS
{
  /**
      @brief Meta information about chemical modification of a sample

      Representation of some kind of modification.
      It hold information about what amino acids are modified and how much the mass changes.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI Modification :
    public SampleTreatment
  {
public:
    /// specificity of the reagent.
    enum SpecificityType
    {
      AA       ///< specified amino acids are modified
      , AA_AT_CTERM      ///< specified amino acids are modified, if they are at the C-terminus
      , AA_AT_NTERM      ///< specified amino acids are modified, if they are at the N-terminus
      , CTERM      ///< the C-terminus is modified
      , NTERM      ///< the N-terminus is modified
      , SIZE_OF_SPECIFICITYTYPE
    };
    /// Names of specificity types
    static const std::string NamesOfSpecificityType[SIZE_OF_SPECIFICITYTYPE];

    /// Default constructor
    Modification();
    /// Copy constructor
    Modification(const Modification &) = default;
    /// Move constructor
    Modification(Modification&&) = default;
    /// Destructor
    ~Modification() override;

    /// Assignment operator
    Modification & operator=(const Modification &) = default;
    /// Move assignment operator
    Modification& operator=(Modification&&) & = default;

    /**
        @brief Equality operator

    Although this operator takes a reference to a SampleTreatment as argument
    it tests for the equality of Tagging instances!
  */
    bool operator==(const SampleTreatment & rhs) const override;

    /// clone method. See SampleTreatment
    SampleTreatment * clone() const override;

    /// returns the name of the reagent that was used (default: "")
    const String & getReagentName() const;
    /// sets the name of the reagent that was used
    void setReagentName(const String & reagent_name);

    /// returns the mass change (default: 0.0)
    double getMass() const;
    /// sets the mass change
    void setMass(double mass);

    /// returns the specificity of the reagent (default: AA)
    const SpecificityType & getSpecificityType() const;
    /// sets the specificity of the reagent
    void setSpecificityType(const SpecificityType & specificity_type);

    /// returns a string containing the one letter code of the amino acids that are affected by the reagent. (default: "")
    const String & getAffectedAminoAcids() const;
    /// returns a string containing the one letter code of the amino acids that are affected by the reagent. Do not separate them by space, tab or comma!
    void setAffectedAminoAcids(const String & affected_amino_acids);

protected:
    String reagent_name_;
    double mass_;
    SpecificityType specificity_type_;
    String affected_amino_acids_;
  };
} // namespace OpenMS

