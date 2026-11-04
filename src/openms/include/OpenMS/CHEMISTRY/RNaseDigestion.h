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

#include <set>

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
    using ConstRibonucleotidePtr = const Ribonucleotide*;

    /// Detailed digestion product including sequence and parent coordinates
    struct DigestionProduct
    {
      NASequence fragment;
      std::pair<Size, Size> position;
    };

    /// Cleavage-sensitive modification groups split by cleavage direction
    struct CleavageSensitiveModGroups
    {
      std::set<ConstRibonucleotidePtr> cuts_before_sensitive;
      std::set<ConstRibonucleotidePtr> cuts_after_sensitive;

      /// Return the union of modifications sensitive to either cleavage direction.
      std::set<ConstRibonucleotidePtr> combined() const
      {
        std::set<ConstRibonucleotidePtr> all = cuts_before_sensitive;
        all.insert(cuts_after_sensitive.begin(), cuts_after_sensitive.end());
        return all;
      }
    };


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

    /**
       @brief Performs the enzymatic digestion of a RNA and returns fragments with parent coordinates

       Only fragments of appropriate length (between @p min_length and @p max_length) are returned.
       Enzyme-specific terminal gains are applied to the reported fragment sequences.

       @param[in] rna RNA parent sequence
       @param[out] output Digestion products with zero-based (start, length) coordinates
       @param[in] min_length Minimum product length
       @param[in] max_length Maximum product length (zero means unlimited)
    */
    void digest(const NASequence& rna, std::vector<DigestionProduct>& output,
                Size min_length = 0, Size max_length = 0) const;

    /**
       @brief Returns the positions of digestion products in the RNA as pairs: (start, length)

       This is useful when callers need to associate digested fragments with parent coordinates.

       @param[in] rna RNA parent sequence
       @param[in] min_length Minimum product length
       @param[in] max_length Maximum product length (zero means unlimited)
    */
    std::vector<std::pair<Size, Size>> getFragmentPositions(
      const NASequence& rna, Size min_length = 0, Size max_length = 0) const;

    /**
       @brief Infer which variable modifications can block cleavage for the configured enzyme

       A modification is classified as cleavage-sensitive if its origin residue matches the
       enzyme cleavage regex at a boundary position, but the modified residue code no longer matches.

       @param[in] variable_modifications Candidate variable modifications
    */
    CleavageSensitiveModGroups inferCleavageSensitiveMods(
      const std::set<ConstRibonucleotidePtr>& variable_modifications) const;

    /**
       @brief Digest RNA while allowing cleavage-sensitive modifications to block adjacent cuts

       Enumerate eligible modification placements in candidate fragments, including internal
       missed-cleavage sites and parent termini. Cleavage boundaries and the number of remaining
       missed cleavages are evaluated on the modified sequence. The number of modifications
       is limited by @p max_sensitive_mods_per_fragment. Enzyme terminal gains are applied
       to all returned fragment sequences.

       @param[in] rna RNA parent sequence
       @param[in] cleavage_sensitive_mods Modifications that block enzyme cuts
       @param[in] max_sensitive_mods_per_fragment Maximum variable modification count
       @param[out] output Digestion products with parent coordinates
       @param[in] min_length Minimum product length
       @param[in] max_length Maximum product length (zero means unlimited)
    */
    void digestWithCleavageSensitiveMods(
      const NASequence& rna,
      const CleavageSensitiveModGroups& cleavage_sensitive_mods,
      Size max_sensitive_mods_per_fragment,
      std::vector<DigestionProduct>& output,
      Size min_length = 0,
      Size max_length = 0) const;

    /// A deduplicated digestion candidate with its evidence in the database sequences of the run.
    struct OPENMS_DLLAPI DigestedOligo
    {
      NASequence sequence;
      std::vector<IdentificationData::SequenceEvidence> sequence_evidence;
      IdentificationData::TargetDecoy target_decoy = IdentificationData::TargetDecoy::UNKNOWN;

      /**
         @brief Add an occurrence of the candidate in a database sequence

         The first occurrence sets the target/decoy state; occurrences in target and decoy sequences make it both,
         and one in a sequence of unknown state makes it unknown. An occurrence that is already known is not added again.

         @param[in] entry The database sequence
         @param[in] rna The sequence of @p entry (its neighbors of the occurrence are taken from it)
         @param[in] position Zero-based start and length of the occurrence in @p rna
      */
      void addEvidence(const IdentificationData::DatabaseSequence& entry, const NASequence& rna, const std::pair<Size, Size>& position);
    };
    /// Digest the database sequences of an oligonucleotide run without creating experimental spectrum matches.
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

    /// Apply enzyme-specific 5'/3' terminal gains to a fragment based on its parent coordinates
    void applyTerminalGains_(NASequence& fragment,
                         const std::pair<Size, Size>& pos,
                         Size parent_size) const;
  };

} // namespace OpenMS
