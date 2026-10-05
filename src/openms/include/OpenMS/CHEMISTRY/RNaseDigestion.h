// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Samuel Wein $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/DATASTRUCTURES/RegularExpression.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

namespace OpenMS
{
  /**
     @brief Class for the enzymatic digestion of RNAs

     @see @ref DigestionEnzymeRNA

     @ingroup Chemistry
  */
  class OPENMS_DLLAPI RNaseDigestion: public EnzymaticDigestion
  {
  public:
    /// Sets the enzyme for the digestion
    void setEnzyme(const DigestionEnzyme* enzyme) override;

    /// Sets the enzyme for the digestion (by name)
    void setEnzyme(const std::string& name);

    /**
       @brief Performs the enzymatic digestion of a (potentially modified) RNA

       Only fragments of appropriate length (between @p min_length and @p max_length) are returned.
    */
    void digest(const NASequence& rna, std::vector<NASequence>& output,
                Size min_length = 0, Size max_length = 0) const;

    /// A deduplicated digestion candidate with owning sequence-to-parent evidence.
    struct DigestedOligo
    {
      NASequence sequence;
      std::vector<IdentificationData::ParentEvidence> parent_evidence;
      IdentificationData::TargetDecoy target_decoy = IdentificationData::TargetDecoy::UNKNOWN;
    };
    /// Digest RNA parents without creating experimental spectrum matches.
    std::vector<DigestedOligo> digest(const IdentificationData::Run& run, Size min_length = 0, Size max_length = 0) const;

  protected:
    const Ribonucleotide* five_prime_gain_; ///< 5' mod added by the enzyme
    const Ribonucleotide* three_prime_gain_; ///< 3' mod added by the enzyme
    std::vector<RegularExpression>
      cuts_after_regexes_; ///< a vector of reg. exp. for enzyme cutting pattern, each regex represents a single nucleotide
    std::vector<RegularExpression> cuts_before_regexes_; ///< a vector reg. exp. for enzyme cutting pattern

    /**
       @brief Returns the positions of digestion products in the RNA as pairs: (start, length)
     */
    std::vector<std::pair<Size, Size>> getFragmentPositions_(
      const NASequence& rna, Size min_length, Size max_length)
      const;
  };

} // namespace OpenMS

