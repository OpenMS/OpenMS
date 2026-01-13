// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ProFormaData.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <iosfwd>
#include <vector>

namespace OpenMS
{

  /**
    @brief Writer for ProForma v2 peptidoform notation

    This class provides serialization of ProForma AST structures back to
    ProForma strings. It supports lossless roundtrip conversion when the
    AST was populated from parsing, preserving the original notation style
    (e.g., using "Oxidation" vs "UNIMOD:35" based on what was originally parsed).

    Output Format Rules:
    - Global mods: `<[MOD]@A,C>` or `<13C>` for isotope replacement
    - N-terminal: `[mod]-SEQUENCE`
    - C-terminal: `SEQUENCE-[mod]`
    - Modifications: `A[mod1][mod2]`
    - Alternatives: `[mod1|mod2]`
    - Cross-links: `[mod#XL1]`
    - Glycan: `[Glycan:HexNAc2Hex5]`
    - CV accession: `[UNIMOD:35]`
    - Mass delta: `[+15.9949]`
    - Formula: `[Formula:C2H3NO]`
    - Charge: `/2` or `/[M+2H]2+`
    - Unlocalised: `[Phospho]?` or `[Phospho]^2?`
    - Labile: `{Glycan:Hex}`
    - Ambiguous region: `(?DQ)`
    - Modified range: `(EOSFORMS)[+19.0523]`
    - Multiple chains: `PEPTIDE//PEPTIDE`

    Usage example:
    @code
    Peptidoform peptidoform;
    // ... populate peptidoform ...
    String proforma_str = ProFormaWriter::toString(peptidoform);

    PeptidoformIon ion;
    // ... populate ion with charge ...
    String ion_str = ProFormaWriter::toString(ion);
    @endcode

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ProFormaWriter
  {
  public:

    /**
      @brief Convert a Peptidoform to a ProForma string

      Serializes a single peptide chain including all modifications,
      but without charge state information.

      @param peptidoform The peptidoform to serialize
      @return The ProForma string representation
    */
    static String toString(const Peptidoform& peptidoform);

    /**
      @brief Convert a PeptidoformIon to a ProForma string

      Serializes one or more peptide chains with optional charge state.
      Multiple chains are separated by "//".

      @param ion The peptidoform ion to serialize
      @return The ProForma string representation including charge if present
    */
    static String toString(const PeptidoformIon& ion);

  private:

    /**
      @brief Write global modification entries to stream

      Handles both IsotopeReplacement (<13C>) and GlobalModification (<[TMT6plex]@K>).

      @param os Output stream
      @param mods Vector of global modification entries
    */
    static void writeGlobalMods_(std::ostream& os, const std::vector<GlobalModEntry>& mods);

    /**
      @brief Write an isotope replacement to stream

      Formats as <ISOTOPE>, e.g., <13C>, <15N>, <D>

      @param os Output stream
      @param isotope The isotope replacement
    */
    static void writeIsotopeReplacement_(std::ostream& os, const IsotopeReplacement& isotope);

    /**
      @brief Write a global modification to stream

      Formats as <[MOD]@TARGETS>, e.g., <[TMT6plex]@K,N-term>

      @param os Output stream
      @param mod The global modification
    */
    static void writeGlobalModification_(std::ostream& os, const GlobalModification& mod);

    /**
      @brief Write a modification (with alternatives) to stream

      Handles single modifications and alternatives separated by |.
      Includes labels like #XL1 or #g1(0.90) if present.

      @param os Output stream
      @param mod The modification with its alternatives
    */
    static void writeModification_(std::ostream& os, const Modification& mod);

    /**
      @brief Write a modification tag (variant) to stream

      Dispatches to the appropriate writer based on tag type:
      CvAccession, NamedMod, MassDelta, FormulaTag, GlycanComposition, InfoTag

      @param os Output stream
      @param tag The modification tag variant
    */
    static void writeModificationTag_(std::ostream& os, const ModificationTag& tag);

    /**
      @brief Write a CV accession to stream

      Formats as DATABASE:ACCESSION, e.g., UNIMOD:35, MOD:00046

      @param os Output stream
      @param cv The CV accession
    */
    static void writeCvAccession_(std::ostream& os, const CvAccession& cv);

    /**
      @brief Write a named modification to stream

      Formats as NAME or PREFIX:NAME, e.g., Oxidation, U:Oxidation

      @param os Output stream
      @param named The named modification
    */
    static void writeNamedMod_(std::ostream& os, const NamedMod& named);

    /**
      @brief Write a mass delta to stream

      Uses original_text if available for lossless roundtrip.
      Otherwise formats as +/-MASS, e.g., +15.9949, -17.03

      @param os Output stream
      @param delta The mass delta
    */
    static void writeMassDelta_(std::ostream& os, const MassDelta& delta);

    /**
      @brief Write a formula tag to stream

      Formats as Formula:FORMULA or Formula:FORMULA:z+CHARGE

      @param os Output stream
      @param formula The formula tag
    */
    static void writeFormulaTag_(std::ostream& os, const FormulaTag& formula);

    /**
      @brief Write a glycan composition to stream

      Formats as Glycan:COMPOSITION, e.g., Glycan:HexNAc2Hex5

      @param os Output stream
      @param glycan The glycan composition
    */
    static void writeGlycanComposition_(std::ostream& os, const GlycanComposition& glycan);

    /**
      @brief Write an info tag to stream

      Formats as INFO:TEXT

      @param os Output stream
      @param info The info tag
    */
    static void writeInfoTag_(std::ostream& os, const InfoTag& info);

    /**
      @brief Write a label to stream

      Formats as #ID or #ID(SCORE), e.g., #XL1, #g1(0.90), #BRANCH

      @param os Output stream
      @param label The label
    */
    static void writeLabel_(std::ostream& os, const Label& label);

    /**
      @brief Write unlocalised modifications to stream

      Formats as [MOD]? or [MOD]^N?, e.g., [Phospho]?, [Phospho]^2?

      @param os Output stream
      @param mods Vector of unlocalised modifications
    */
    static void writeUnlocalisedMods_(std::ostream& os, const std::vector<UnlocalisedMod>& mods);

    /**
      @brief Write labile modifications to stream

      Formats as {MOD}, e.g., {Glycan:Hex}

      @param os Output stream
      @param mods Vector of labile modifications
    */
    static void writeLabileModifications_(std::ostream& os, const std::vector<LabileModification>& mods);

    /**
      @brief Write N-terminal modifications to stream

      Formats as [MOD]-, e.g., [Acetyl]-

      @param os Output stream
      @param mods Vector of N-terminal modifications
    */
    static void writeNTermMods_(std::ostream& os, const std::vector<Modification>& mods);

    /**
      @brief Write C-terminal modifications to stream

      Formats as -[MOD], e.g., -[Amidated]

      @param os Output stream
      @param mods Vector of C-terminal modifications
    */
    static void writeCTermMods_(std::ostream& os, const std::vector<Modification>& mods);

    /**
      @brief Write the sequence with modifications to stream

      Handles SequenceElement, AmbiguousRegion, and ModifiedRange.

      @param os Output stream
      @param seq Vector of sequence sections
    */
    static void writeSequence_(std::ostream& os, const std::vector<SequenceSection>& seq);

    /**
      @brief Write a single sequence element to stream

      Formats as AA[MOD][MOD]..., e.g., M[Oxidation], K[TMT6plex]

      @param os Output stream
      @param elem The sequence element
    */
    static void writeSequenceElement_(std::ostream& os, const SequenceElement& elem);

    /**
      @brief Write an ambiguous region to stream

      Formats as (?AA...), e.g., (?DQ)

      @param os Output stream
      @param region The ambiguous region
    */
    static void writeAmbiguousRegion_(std::ostream& os, const AmbiguousRegion& region);

    /**
      @brief Write a modified range to stream

      Formats as (SEQ)[MOD], e.g., (EOSFORMS)[+19.0523]

      @param os Output stream
      @param range The modified range
    */
    static void writeModifiedRange_(std::ostream& os, const ModifiedRange& range);

    /**
      @brief Write a charge state to stream

      Formats as /CHARGE or /[ADDUCTS]CHARGE+, e.g., /2, /[M+2H]2+

      @param os Output stream
      @param charge The charge state
    */
    static void writeChargeState_(std::ostream& os, const ChargeState& charge);

    /**
      @brief Get the string representation of a CV database prefix

      @param db The CV database enum
      @return The prefix string (e.g., "UNIMOD", "MOD", "XLMOD")
    */
    static const char* cvDatabaseToString_(CvDatabase db);

    /**
      @brief Get the single-letter hint for a CV database

      @param db The CV database enum
      @return The single-letter hint (e.g., "U", "M", "X")
    */
    static const char* cvDatabaseToHint_(CvDatabase db);

    /**
      @brief Get the source prefix string for mass delta

      @param source The mass delta source enum
      @return The source prefix (e.g., "Obs", "U", "M") or empty for NONE
    */
    static const char* massSourceToString_(MassDelta::Source source);
  };

} // namespace OpenMS
