// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

