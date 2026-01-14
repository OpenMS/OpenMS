// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ProFormaData.h>
#include <OpenMS/CHEMISTRY/ProFormaError.h>
#include <OpenMS/CHEMISTRY/ProFormaTokenizer.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <string>
#include <string_view>
#include <vector>

namespace OpenMS
{

  // Forward declaration
  class AASequence;

  /**
    @brief Recursive descent parser for ProForma v2 peptidoform notation

    This class parses ProForma strings into an Abstract Syntax Tree (AST) representation.
    The AST structures are defined in ProFormaData.h.

    The parser implements the ProForma v2 grammar:
    - peptidoform_ion -> peptidoform [charge_state | adduct_ion]
    - peptidoform -> [global_mods] [n_term_mod "-"] sequence ["-" c_term_mod]
    - sequence -> {amino_acid [modification_list]}+
    - modification_list -> "[" modification {"," modification}* "]"
    - modification -> named_mod | mass_delta | formula | cv_accession | glycan | info_tag
    - global_mods -> "<" global_mod {"," global_mod}* ">"
    - cross_link_label -> "#" identifier ["(" occurrence ")"]

    Usage example:
    @code
    String input = "EM[UNIMOD:35]K";
    Peptidoform pf = ProFormaParser::parse(input);
    // pf now contains the parsed AST
    @endcode

    @see ProFormaData.h for AST structure definitions
    @see ProFormaTokenizer for lexical analysis
    @see ProFormaError for error handling

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ProFormaParser
  {
  public:

    /**
      @brief Parse a ProForma string into a Peptidoform AST

      This is the main entry point for parsing simple peptidoforms without
      charge state information.

      @param input The ProForma string to parse
      @return The parsed Peptidoform AST
      @throws ProFormaParseError if the input is invalid

      @note For peptidoforms with charge state (e.g., "PEPTIDE/2"), use parseIon()
    */
    static Peptidoform parse(const String& input);

    /**
      @brief Parse a ProForma string into a PeptidoformIon AST

      This entry point handles the full ProForma notation including:
      - Multiple chains (separated by //)
      - Charge state specification (/2, /+2, /-1)
      - Adduct ion specification (/[Na:z+1,H:z+1])

      @param input The ProForma string to parse
      @return The parsed PeptidoformIon AST
      @throws ProFormaParseError if the input is invalid
    */
    static PeptidoformIon parseIon(const String& input);

    /**
      @brief Convert a Peptidoform AST back to ProForma string notation

      @param pf The Peptidoform to convert
      @param mode Write mode: LOSSLESS preserves original formatting, CANONICAL produces normalized output
      @return The ProForma string representation
    */
    static String toString(const Peptidoform& pf,
                           ProFormaWriteMode mode = ProFormaWriteMode::LOSSLESS);

    /**
      @brief Convert a PeptidoformIon AST back to ProForma string notation

      @param pfi The PeptidoformIon to convert
      @param mode Write mode: LOSSLESS preserves original formatting, CANONICAL produces normalized output
      @return The ProForma string representation
    */
    static String toString(const PeptidoformIon& pfi,
                           ProFormaWriteMode mode = ProFormaWriteMode::LOSSLESS);

    // ---- AASequence Conversion Methods ----

    /**
      @brief Resolve all modifications in a Peptidoform using ModificationsDB

      Looks up each modification tag (CV accession, named mod, mass delta) in
      ModificationsDB and stores the resolved ResidueModification pointer.

      @param pf The Peptidoform to resolve (modified in place)
      @note Modifications that cannot be resolved will have resolved_mod = nullptr
    */
    static void resolveModifications(Peptidoform& pf);

    /**
      @brief Convert a Peptidoform to an OpenMS AASequence

      @param pf The Peptidoform to convert
      @param policy How to handle unconvertible modifications
      @return The equivalent AASequence
      @throws Exception::ConversionError if STRICT policy and conversion not possible

      @note Call resolveModifications() first, or this will resolve automatically
    */
    static AASequence toAASequence(
      const Peptidoform& pf,
      AASequenceConversionPolicy policy = AASequenceConversionPolicy::STRICT);

    /**
      @brief Create a Peptidoform from an OpenMS AASequence

      Converts an AASequence with modifications to ProForma notation.
      Uses CV accessions (UNIMOD) where available, otherwise named modifications.

      @param seq The AASequence to convert
      @return The equivalent Peptidoform AST
    */
    static Peptidoform fromAASequence(const AASequence& seq);

    /**
      @brief Check if a Peptidoform can be fully represented as an AASequence

      Returns true if all modifications can be resolved and there are no
      unsupported features (ambiguous regions, cross-links, etc.)

      @param pf The Peptidoform to check
      @return True if conversion is possible without issues
    */
    static bool isRepresentableAsAASequence(const Peptidoform& pf);

    /**
      @brief Get a list of all issues that would arise during AASequence conversion

      Returns detailed information about every aspect of the Peptidoform that
      cannot be represented in an AASequence.

      @param pf The Peptidoform to analyze
      @return Vector of conversion issues (empty if fully convertible)
    */
    static std::vector<ConversionIssue> getAASequenceConversionIssues(const Peptidoform& pf);

  private:

    /// Private constructor - use static methods
    explicit ProFormaParser(std::string_view input);

    // ---- High-level parsing methods ----

    /// Parse a complete PeptidoformIon (multiple chains + charge)
    PeptidoformIon parsePeptidoformIon_();

    /// Parse a single Peptidoform (one chain)
    Peptidoform parsePeptidoform_();

    /// Parse global modifications: < ... >
    std::vector<GlobalModEntry> parseGlobalMods_();

    /// Parse a single global modification entry
    GlobalModEntry parseGlobalModEntry_();

    /// Parse isotope replacement: <13C>, <15N>, <D>
    IsotopeReplacement parseIsotopeReplacement_();

    /// Parse global modification with locations: <[mod]@locations>
    GlobalModification parseGlobalModification_();

    /// Parse unlocalised modifications: [mod]?
    std::vector<UnlocalisedMod> parseUnlocalisedMods_();

    /// Parse labile modifications: {mod}
    std::vector<LabileModification> parseLabileModifications_();

    /// Parse the amino acid sequence with modifications
    std::vector<SequenceSection> parseSequence_();

    /// Parse a single sequence element (amino acid + mods)
    SequenceElement parseSequenceElement_();

    /// Parse an ambiguous region: (?XY)
    AmbiguousRegion parseAmbiguousRegion_();

    /// Parse a modified range: (XYZ)[mod]
    ModifiedRange parseModifiedRange_();

    /// Parse terminal modifications: [mod1][mod2]...
    std::vector<Modification> parseTerminalMods_();

    // ---- Modification parsing ----

    /// Parse a modification list: [mod1, mod2, ...]
    std::vector<Modification> parseModificationList_();

    /// Parse a single modification (may have alternatives with |)
    Modification parseModification_();

    /// Parse a single modification tag (no alternatives)
    std::pair<ModificationTag, std::optional<Label>> parseModificationTagWithLabel_();

    /// Parse a modification tag
    ModificationTag parseModificationTag_();

    /// Parse a named modification: Oxidation, U:Oxidation
    NamedMod parseNamedMod_();

    /// Parse a named modification with a known CV hint prefix
    NamedMod parseNamedMod_(char cv_hint);

    /// Parse a CV accession: UNIMOD:35, MOD:00046
    CvAccession parseCvAccession_();

    /// Parse a mass delta: +15.9949, Obs:+79.978
    MassDelta parseMassDelta_();

    /// Parse a formula tag: Formula:C12H20O2
    FormulaTag parseFormulaTag_();

    /// Parse a glycan composition: Glycan:HexNAc1Hex2
    GlycanComposition parseGlycanComposition_();

    /// Parse an info tag: INFO:text
    InfoTag parseInfoTag_();

    /// Parse a label: #XL1, #BRANCH, #g1(0.90)
    Label parseLabel_();

    // ---- Charge state parsing ----

    /// Parse charge state: /2, /+2, /[Na:z+1]
    std::optional<ChargeState> parseChargeState_();

    /// Parse adduct ions: [Na:z+1, H:z+1]
    std::vector<AdductIon> parseAdductIons_();

    /// Parse a single adduct ion: Na:z+1
    AdductIon parseAdductIon_();

    // ---- Helper methods ----

    /// Get the current token
    ProFormaTokenizer::Token current_();

    /// Look at the next token without consuming
    ProFormaTokenizer::Token peek_();

    /// Consume and return the current token
    ProFormaTokenizer::Token advance_();

    /// Check if current token matches expected type
    bool check_(ProFormaTokenizer::TokenType type);

    /// Check if current token matches expected type, consume if true
    bool match_(ProFormaTokenizer::TokenType type);

    /// Expect a specific token type, throw error if not found
    ProFormaTokenizer::Token expect_(ProFormaTokenizer::TokenType type, const char* expected_desc);

    /// Check if at end of input
    bool isAtEnd_();

    /// Throw a parse error at the current position
    [[noreturn]] void error_(ProFormaErrorCode code, const char* message);

    /// Throw a parse error at a specific position
    [[noreturn]] void errorAt_(ProFormaErrorCode code, size_t pos, const char* message);

    /// Parse a CV database prefix from identifier
    std::optional<CvDatabase> parseCvDatabasePrefix_(const std::string_view& id);

    /// Check if identifier is a valid amino acid
    static bool isAminoAcid_(char c);

    /// Check if the current position could start a modification tag content
    bool looksLikeModificationTagContent_();

    /// Check if current position has N-terminal modification pattern ([mod]-)
    bool hasNTerminalModPattern_();

    // ---- Member variables ----

    /// The tokenizer for lexical analysis
    ProFormaTokenizer tokenizer_;

    /// The original input string (for error messages)
    std::string input_;

    /// Current token (cached)
    ProFormaTokenizer::Token current_token_;

    /// Whether we have a cached current token
    bool has_current_ = false;
  };

} // namespace OpenMS
