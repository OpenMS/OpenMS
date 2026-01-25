// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace OpenMS
{

  // Forward declarations
  class ResidueModification;
  class AASequence;
  class MSSpectrum;

  /**
    @brief ProForma v2 peptidoform notation parser and data structures

    This class provides parsing, serialization, and conversion functionality for
    the ProForma v2 peptidoform notation standard. It contains nested types that
    form the Abstract Syntax Tree (AST) representation of parsed ProForma strings.

    Usage example:
    @code
    String input = "EM[UNIMOD:35]K";
    ProForma::Peptidoform pf = ProForma::parse(input);
    // pf now contains the parsed AST

    // Convert back to string
    String output = ProForma::toString(pf);

    // Convert to AASequence
    AASequence seq = ProForma::toAASequence(pf);
    @endcode

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI ProForma
  {
  public:
    //==========================================================================
    // Enums
    //==========================================================================

    /**
      @brief Conversion policy for transforming Peptidoform to AASequence

      Controls how the conversion handles modifications that cannot be directly
      represented in AASequence (e.g., unlocalised, labile, or ambiguous modifications).
    */
    enum class ConversionPolicy
    {
      FAIL_ON_LOSS,      ///< Fail if any modification cannot be fully represented
      DROP_UNLOCALISED,  ///< Drop unlocalised, labile, and global modifications
      BEST_EFFORT        ///< Try to convert as much as possible, skip unsupported
    };


    /**
      @brief Issue type for AASequence conversion problems
    */
    enum class ConversionIssueType
    {
      UNRESOLVED_MOD,      ///< Modification could not be found in ModificationsDB
      UNLOCALISED_MOD,     ///< Modification has no specific position
      LABILE_MOD,          ///< Labile modification (lost during fragmentation)
      GLOBAL_MOD,          ///< Global modification (applies to multiple sites)
      AMBIGUOUS_MOD,       ///< Ambiguously localized modification
      AMBIGUOUS_REGION,    ///< Ambiguous amino acid region
      MODIFIED_RANGE,      ///< Modified range (position uncertain)
      CROSS_LINK,          ///< Cross-link between chains
      MULTIPLE_CHAINS,     ///< Multiple peptide chains
      ALTERNATIVE_MODS,    ///< Multiple alternative modifications (|)
      UNSUPPORTED_FEATURE  ///< Other unsupported ProForma feature
    };


    /**
      @brief Write mode for ProForma string serialization

      Controls whether the output preserves original formatting (LOSSLESS) or
      produces a normalized, deterministic output (CANONICAL).
    */
    enum class WriteMode
    {
      LOSSLESS,   ///< Preserve original spelling/formatting where possible (e.g., mass delta text)
      CANONICAL   ///< Normalized output: uppercase CV prefixes, sorted mods, 4 decimal places for masses
    };


    /**
      @brief Controlled vocabulary database prefix for modification accessions

      Identifies the source database for a modification accession in ProForma notation.
      Examples: UNIMOD:35, MOD:00046, XLMOD:02001, GNO:G59626AS
    */
    enum class CvDatabase
    {
      UNIMOD, ///< UniMod database (https://www.unimod.org/)
      MOD,    ///< PSI-MOD ontology (https://www.ebi.ac.uk/ols/ontologies/mod)
      RESID,  ///< RESID database
      XLMOD,  ///< Cross-linking modifications ontology
      GNO     ///< Glycan naming ontology
    };


    /**
      @brief Error codes for programmatic handling of ProForma parse errors.

      These error codes provide machine-readable categorization of parsing failures,
      enabling downstream code to handle specific error types appropriately.
    */
    enum class ErrorCode
    {
      UNEXPECTED_CHARACTER,       ///< Unexpected character encountered during parsing
      UNCLOSED_BRACKET,           ///< Opening bracket without matching close bracket
      UNMATCHED_BRACKET,          ///< Closing bracket without matching open bracket
      INVALID_CV_PREFIX,          ///< Invalid controlled vocabulary prefix (e.g., not UNIMOD, MOD, etc.)
      INVALID_CV_ACCESSION,       ///< Invalid CV accession number format
      INVALID_AMINO_ACID,         ///< Invalid amino acid one-letter code
      INVALID_MASS_VALUE,         ///< Invalid mass value format or value
      INVALID_FORMULA,            ///< Invalid chemical formula
      UNKNOWN_MONOSACCHARIDE,     ///< Unknown monosaccharide abbreviation
      DANGLING_CROSSLINK_LABEL,   ///< Crosslink label without a matching partner
      EMPTY_SEQUENCE,             ///< Empty sequence provided
      INVALID_CHARGE,             ///< Invalid charge state specification
      INVALID_OCCURRENCE_SPECIFIER, ///< Invalid occurrence specifier (e.g., ^2)
      UNEXPECTED_END_OF_INPUT,    ///< Unexpected end of input string
      INTERNAL_ERROR              ///< Internal parser error (should not occur)
    };


    //==========================================================================
    // Data Structures
    //==========================================================================

    /**
      @brief Description of a conversion issue from Peptidoform to AASequence

      Records problems encountered when attempting to convert a ProForma
      Peptidoform to an OpenMS AASequence representation.
    */
    struct OPENMS_DLLAPI ConversionIssue
    {
      ConversionIssueType type;   ///< The type of issue
      String description;         ///< Human-readable description
      size_t position;            ///< Position in sequence (SIZE_MAX if not position-specific)
    };


    /**
      @brief Controlled vocabulary accession for a modification

      Represents a modification specified by a CV accession number, e.g., UNIMOD:35 for Oxidation.
      The accession string contains only the identifier portion (e.g., "35" for UNIMOD:35).
    */
    struct OPENMS_DLLAPI CvAccession
    {
      CvDatabase database;  ///< The source database (UNIMOD, MOD, RESID, XLMOD, or GNO)
      String accession;     ///< The accession identifier (e.g., "35" for UNIMOD:35, full string for GNO)
    };


    /**
      @brief Named modification with optional CV prefix hint

      Represents a modification specified by name, optionally with a CV prefix hint
      to disambiguate which database to search (e.g., "U:Oxidation" for UniMod,
      "M:Oxidation" for PSI-MOD).
    */
    struct OPENMS_DLLAPI NamedMod
    {
      std::optional<CvDatabase> cv_hint;  ///< Optional CV prefix hint (U, M, R, X, G)
      String name;                        ///< The modification name (e.g., "Oxidation", "Phospho")
    };


    /**
      @brief Mass delta modification with optional source hint

      Represents a modification specified by mass difference. The source hint indicates
      which database was used to observe/calculate this mass (e.g., Obs:+79.978).
      The original_text preserves the exact formatting for lossless roundtrip.
    */
    struct OPENMS_DLLAPI MassDelta
    {
      /// Source hint for mass delta values
      enum class Source
      {
        NONE,  ///< No source hint specified
        OBS,   ///< Observed mass (Obs:)
        U,     ///< UniMod lookup (U:)
        M,     ///< PSI-MOD lookup (M:)
        R,     ///< RESID lookup (R:)
        X,     ///< XLMOD lookup (X:)
        G      ///< GNO lookup (G:)
      };

      Source source = Source::NONE;  ///< Optional source hint prefix
      double mass;                   ///< The mass delta value in Daltons
      String original_text;          ///< Original text for lossless roundtrip (e.g., "+15.99" vs "+15.9900")
    };


    /**
      @brief Chemical formula with optional charge

      Represents a modification specified by chemical formula. The optional charge
      is specified via the :z+N suffix in ProForma (e.g., Formula:C12H20O2:z+2).
    */
    struct OPENMS_DLLAPI FormulaTag
    {
      String formula_string;        ///< The chemical formula string (e.g., "C12H20O2")
      std::optional<int> charge;    ///< Optional charge from :z+N suffix
    };


    /**
      @brief Glycan composition specification

      Represents a glycan modification as a composition of monosaccharides.
      Each component can be either a named monosaccharide (e.g., "Hex", "HexNAc")
      or a custom formula specification.

      Example: Glycan:HexNAc1Hex2 -> [(HexNAc, 1), (Hex, 2)]
    */
    struct OPENMS_DLLAPI GlycanComposition
    {
      /// A monosaccharide component: either a name (String) or a custom formula (FormulaTag)
      using Monosaccharide = std::variant<String, FormulaTag>;

      std::vector<std::pair<Monosaccharide, int>> components;  ///< List of (monosaccharide, count) pairs
    };


    /**
      @brief Info tag for arbitrary text annotations

      Represents an INFO: tag in ProForma that carries arbitrary text metadata
      about a modification or site. Example: INFO:provenance_data
    */
    struct OPENMS_DLLAPI InfoTag
    {
      String text;  ///< The info text content
    };


    /**
      @brief Position constraint specifying allowed residues for a modification

      Represents a Position: tag in ProForma that constrains where a modification
      can be localized. This is typically used as an alternative to a modification
      to indicate its possible sites.

      Example: [Oxidation|Position:M] means Oxidation can only occur at M residues
    */
    struct OPENMS_DLLAPI PositionConstraint
    {
      std::vector<char> residues;  ///< List of allowed amino acid residues
      bool n_term = false;         ///< True if modification can be at N-terminus
      bool c_term = false;         ///< True if modification can be at C-terminus
    };


    /**
      @brief Variant type representing any modification tag content

      A ModificationTag can be one of:
      - CvAccession: CV database accession (e.g., UNIMOD:35)
      - NamedMod: Named modification with optional CV hint (e.g., Oxidation, U:Oxidation)
      - MassDelta: Mass difference with optional source (e.g., +15.9949, Obs:+79.978)
      - FormulaTag: Chemical formula (e.g., Formula:C12H20O2)
      - GlycanComposition: Glycan composition (e.g., Glycan:HexNAc1Hex2)
      - InfoTag: Info annotation (e.g., INFO:comment)
      - PositionConstraint: Allowed residue positions (e.g., Position:MKC)
    */
    using ModificationTag = std::variant<
      CvAccession,
      NamedMod,
      MassDelta,
      FormulaTag,
      GlycanComposition,
      InfoTag,
      PositionConstraint
    >;


    /**
      @brief Label for cross-links, branches, or ambiguous grouping

      Labels are used to:
      - Connect cross-linked sites: `[XLMOD:02001#XL1]...[#XL1]`
      - Mark branch points: `[#BRANCH]`
      - Group ambiguously localized modifications: `[Phospho#g1(0.90)]`
    */
    struct OPENMS_DLLAPI Label
    {
      /// The type of label
      enum class Type
      {
        CROSSLINK,  ///< Cross-link label (e.g., `#XL1`)
        BRANCH,     ///< Branch point label (`#BRANCH`)
        AMBIGUOUS   ///< Ambiguous localization group (e.g., `#g1`)
      };

      Type type;                        ///< The label type
      String identifier;                ///< The label identifier (e.g., `XL1`, `BRANCH`, `g1`)
      std::optional<double> score;      ///< Optional localization score for ambiguous labels (e.g., 0.90)
    };


    /**
      @brief A modification with one or more alternative tags

      In ProForma, a modification can have multiple alternatives separated by |,
      representing uncertainty about the exact modification. Each alternative
      consists of a tag and an optional label.

      Example: K[Phospho|+79.97] has two alternatives

      The resolved_mod field is populated by resolveModifications() and points
      to the ResidueModification in ModificationsDB (for the first/primary alternative).
    */
    struct OPENMS_DLLAPI Modification
    {
      /// Each alternative is a (tag, optional_label) pair
      std::vector<std::pair<ModificationTag, std::optional<Label>>> alternatives;

      /// Resolved modification pointer (populated by resolveModifications)
      /// Points to the ResidueModification for the first alternative, if found
      const ResidueModification* resolved_mod = nullptr;
    };


    /**
      @brief A single amino acid with its modifications

      Represents one position in the peptide sequence: the amino acid residue
      and zero or more modifications attached to it.
    */
    struct OPENMS_DLLAPI SequenceElement
    {
      char amino_acid;                          ///< Single-letter amino acid code (A-Z)
      std::vector<Modification> modifications;  ///< Modifications at this position
    };


    /**
      @brief Ambiguous amino acid region

      Represents a region where the amino acid sequence is uncertain.
      ProForma notation: (?DQ) means either D or Q at this position.
    */
    struct OPENMS_DLLAPI AmbiguousRegion
    {
      std::vector<SequenceElement> elements;  ///< The ambiguous amino acid possibilities
    };


    /**
      @brief Modified sequence range with shared modifications

      Represents a subsequence where one or more modifications apply to the
      entire range, but the exact position is uncertain.

      ProForma notation: (EOSFORMS)[+19.0523] means +19.0523 applies somewhere in EOSFORMS
    */
    struct OPENMS_DLLAPI ModifiedRange
    {
      std::vector<SequenceElement> elements;    ///< The amino acids in the range
      std::vector<Modification> modifications;  ///< Modifications applying to the entire range
    };


    /**
      @brief Variant type representing a section of the sequence

      A SequenceSection can be:
      - SequenceElement: A single amino acid with its modifications
      - AmbiguousRegion: Uncertain amino acid identity (?XY)
      - ModifiedRange: Subsequence with range-wide modifications (XYZ)[mod]
    */
    using SequenceSection = std::variant<
      SequenceElement,
      AmbiguousRegion,
      ModifiedRange
    >;


    /**
      @brief Unlocalised modification with optional occurrence count

      Represents a modification that is known to exist on the peptide but
      whose exact position is unknown. The occurrence specifies how many
      instances of this modification are present.

      ProForma notation: [Phospho]?PEPTIDE or [Phospho]^2?PEPTIDE
    */
    struct OPENMS_DLLAPI UnlocalisedMod
    {
      std::vector<Modification> modifications;  ///< The unlocalised modification(s)
      std::optional<int> occurrence;            ///< Optional occurrence count from ^N suffix
    };


    /**
      @brief Labile modification that may be lost during fragmentation

      Labile modifications are typically lost during ionization or fragmentation
      and thus may not be observed in MS/MS spectra.

      ProForma notation: {Glycan:Hex}PEPTIDE
    */
    struct OPENMS_DLLAPI LabileModification
    {
      Modification modification;  ///< The labile modification
    };


    /**
      @brief Global modification applied to specific locations

      A global modification applies the same modification to all occurrences
      of specified residues or termini throughout the peptide.

      ProForma notation: `<[TMT6plex]@K,N-term>`
    */
    struct OPENMS_DLLAPI GlobalModification
    {
      Modification modification;        ///< The modification to apply
      std::vector<String> locations;    ///< Target locations ("K", "N-term", "C-term:K", etc.)
    };


    /**
      @brief Isotope replacement for stable isotope labeling

      Represents global replacement of an element with a specific isotope,
      used for stable isotope labeling experiments.

      ProForma notation: `<13C>` or `<15N>` or `<D>`
    */
    struct OPENMS_DLLAPI IsotopeReplacement
    {
      String isotope;  ///< The isotope specification (e.g., "13C", "15N", "D")
    };


    /**
      @brief Variant type for global modification entries

      A GlobalModEntry can be either:
      - IsotopeReplacement: Global isotope substitution (`<13C>`)
      - GlobalModification: Position-specific global mod (`<[TMT6plex]@K>`)
    */
    using GlobalModEntry = std::variant<
      IsotopeReplacement,
      GlobalModification
    >;


    /**
      @brief Adduct ion specification for charge state

      Represents an adduct ion contributing to the charge state of a
      peptidoform ion. Multiple adducts can combine to give the total charge.

      ProForma notation: Na:z+1 in /[Na:z+1,H:z+1]
    */
    struct OPENMS_DLLAPI AdductIon
    {
      String formula;                   ///< The adduct formula (e.g., "Na", "H", "K")
      int charge;                       ///< The charge contribution of this adduct
      std::optional<int> occurrence;    ///< Optional occurrence count from ^N suffix
    };


    /**
      @brief Charge state specification

      The charge state can be specified as either:
      - A simple integer charge (/2, /+2, /-1)
      - A list of adduct ions (/[Na:z+1,H:z+1])
    */
    using ChargeState = std::variant<
      int,                       ///< Simple charge: /2, /+2, /-1
      std::vector<AdductIon>     ///< Adduct list: /[Na:z+1,H:z+1]
    >;


    /**
      @brief A single peptidoform (one peptide chain)

      Represents a complete peptide chain including:
      - Optional name identifier (from v2.1 extension)
      - Global modifications (`<13C>`, `<[TMT6plex]@K>`)
      - Unlocalised modifications ([Phospho]?)
      - Labile modifications ({Glycan:Hex})
      - N-terminal modifications ([Acetyl]-)
      - The amino acid sequence with modifications
      - C-terminal modifications (-[Amidated])
    */
    struct OPENMS_DLLAPI Peptidoform
    {
      std::optional<String> name;                       ///< Optional name from (>name) v2.1 extension
      std::vector<GlobalModEntry> global_mods;          ///< Global modifications: `<13C>`, `<[TMT6plex]@K>`
      std::vector<UnlocalisedMod> unlocalised_mods;     ///< Unlocalised modifications: [Phospho]?
      std::vector<LabileModification> labile_mods;      ///< Labile modifications: {Glycan:Hex}
      std::vector<Modification> n_term_mods;            ///< N-terminal modifications: [Acetyl]-
      std::vector<SequenceSection> sequence;            ///< The sequence with modifications
      std::vector<Modification> c_term_mods;            ///< C-terminal modifications: -[Amidated]
      std::optional<ChargeState> charge;                ///< Optional per-chain charge (for chimeric spectra)
    };


    /**
      @brief A peptidoform ion (one or more chains with optional charge)

      Represents one or more peptide chains that form a single ion species.
      Multiple chains can be present in cross-linked or multi-chain entities.

      ProForma notation: chains are separated by //
    */
    struct OPENMS_DLLAPI PeptidoformIon
    {
      std::optional<String> name;           ///< Optional name from (>>name) v2.1 extension
      std::vector<Peptidoform> chains;      ///< One or more peptide chains (separated by // or + in ProForma)
      std::optional<ChargeState> charge;    ///< Optional charge state specification
      bool is_chimeric = false;             ///< True if chains are chimeric (+ separator), false if cross-linked (//)
    };


    /**
      @brief Cross-link group connecting sites across chains

      Groups together all sites that share a cross-link label. Each site is
      identified by its chain index and position within that chain.

      Derived during parsing from matching `#XL` labels.
    */
    struct OPENMS_DLLAPI CrossLinkGroup
    {
      String label;                                     ///< The cross-link label (e.g., `XL1`)
      std::vector<std::pair<size_t, size_t>> sites;     ///< (chain_index, site_index) pairs
    };


    //==========================================================================
    // ParseError Exception
    //==========================================================================

    /**
      @brief Structured parse error with context for ProForma parsing.

      This exception extends Exception::ParseError to provide more detailed
      information about ProForma parsing failures, including:
      - Machine-readable error code for programmatic handling
      - Position of the error in the input string
      - Context snippets showing surrounding text
      - Expected vs. found information for clearer diagnostics

      Usage example:
      @code
      try {
        ProForma::parse("PEPTIDE[Oxidation");
      } catch (const ProForma::ParseError& e) {
        std::cerr << e.getFormattedMessage() << std::endl;
        // Output:
        // ProForma parse error at position 7: Unclosed bracket
        // Context: PEPTIDE>>>[<<<Oxidation
        // Expected: ']'
        // Found: end of input
      }
      @endcode
    */
    class OPENMS_DLLAPI ParseError :
      public Exception::ParseError
    {
    public:
      /**
        @brief Constructs a ParseError with full context information.

        @param[in] file Source file where the exception was thrown (use __FILE__).
        @param[in] line Source line where the exception was thrown (use __LINE__).
        @param[in] function Function name where the exception was thrown (use OPENMS_PRETTY_FUNCTION).
        @param[in] error_code Machine-readable error code categorizing the error.
        @param[in] error_position Byte position (0-indexed) in the input where the error occurred.
                                  Will be clamped to input.size() if out of bounds.
        @param[in] input The complete input string that was being parsed.
        @param[in] message Human-readable description of the error.
      */
      ParseError(
        const char* file,
        int line,
        const char* function,
        ErrorCode error_code,
        size_t error_position,
        const String& input,
        const String& message
      ) noexcept;

      /// Get the error code for programmatic handling
      ErrorCode getErrorCode() const noexcept
      {
        return code_;
      }

      /// Get byte position of error (0-indexed)
      size_t getPosition() const noexcept
      {
        return position_;
      }

      /// Get context before error (~20 chars)
      const String& getContextBefore() const noexcept
      {
        return context_before_;
      }

      /// Get context after error (~20 chars)
      const String& getContextAfter() const noexcept
      {
        return context_after_;
      }

      /// Get what was expected (if set)
      const String& getExpected() const noexcept
      {
        return expected_;
      }

      /// Get what was found (if set)
      const String& getFound() const noexcept
      {
        return found_;
      }

      /**
        @brief Get formatted error message with context visualization.

        Returns a multi-line string formatted for human reading, showing:
        - Error description with position
        - Context snippet with caret pointing to error location
        - Expected vs. found information (if available)

        @return Formatted error message.
      */
      String getFormattedMessage() const;

      /**
        @brief Set expected/found information for more detailed error messages.

        @param[in] expected Description of what was expected at the error position.
        @param[in] found Description of what was actually found.
      */
      void setExpectedFound(const String& expected, const String& found);

    private:
      ErrorCode code_;              ///< Machine-readable error code
      size_t position_;             ///< Position in input where error occurred
      String context_before_;       ///< Text before the error position
      String context_after_;        ///< Text after the error position
      String expected_;             ///< What was expected (optional)
      String found_;                ///< What was found (optional)

      /**
        @brief Extract context snippets from input around the error position.

        Extracts approximately 20 characters before and after the error position
        for display in error messages. The position is assumed to be already
        clamped to input.size().

        @param[in] input The complete input string.
        @param[in] pos The error position in the input (clamped to input.size()).
      */
      void extractContext_(const String& input, size_t pos);
    };


    //==========================================================================
    // Static Methods - Parsing
    //==========================================================================

    /**
      @brief Parse a ProForma string into a Peptidoform AST

      This is the main entry point for parsing simple peptidoforms without
      charge state information.

      @param[in] input The ProForma string to parse
      @return The parsed Peptidoform AST
      @throws ParseError if the input is invalid

      @note For peptidoforms with charge state (e.g., "PEPTIDE/2"), use parseIon()
    */
    static Peptidoform parse(const String& input);

    /**
      @brief Parse a ProForma string into a PeptidoformIon AST

      This entry point handles the full ProForma notation including:
      - Multiple chains (separated by //)
      - Charge state specification (/2, /+2, /-1)
      - Adduct ion specification (/[Na:z+1,H:z+1])

      @param[in] input The ProForma string to parse
      @return The parsed PeptidoformIon AST
      @throws ParseError if the input is invalid
    */
    static PeptidoformIon parseIon(const String& input);


    //==========================================================================
    // Static Methods - Serialization
    //==========================================================================

    /**
      @brief Convert a Peptidoform AST back to ProForma string notation

      @param[in] pf The Peptidoform to convert
      @param[in] mode Write mode: LOSSLESS preserves original formatting, CANONICAL produces normalized output
      @return The ProForma string representation
    */
    static String toString(const Peptidoform& pf,
                           WriteMode mode = WriteMode::LOSSLESS);

    /**
      @brief Convert a PeptidoformIon AST back to ProForma string notation

      @param[in] pfi The PeptidoformIon to convert
      @param[in] mode Write mode: LOSSLESS preserves original formatting, CANONICAL produces normalized output
      @return The ProForma string representation
    */
    static String toString(const PeptidoformIon& pfi,
                           WriteMode mode = WriteMode::LOSSLESS);


    //==========================================================================
    // Static Methods - JSON Serialization
    //==========================================================================

    /**
      @brief Convert a Peptidoform to JSON string representation

      Serializes the complete Peptidoform AST including all modifications,
      terminal modifications, global modifications, and labels.

      @param[in] pf The Peptidoform to serialize
      @return JSON string representation of the Peptidoform
    */
    static String toJSON(const Peptidoform& pf);

    /**
      @brief Construct a Peptidoform from JSON string

      Deserializes a JSON string back into a Peptidoform AST.

      @param[in] json_str JSON string representation of a Peptidoform
      @return The deserialized Peptidoform
      @throws Exception::ParseError if the JSON is malformed or missing required fields
    */
    static Peptidoform peptidoformFromJSON(const String& json_str);

    /**
      @brief Convert a PeptidoformIon to JSON string representation

      Serializes the complete PeptidoformIon AST including all chains,
      charge state, adduct ions, and cross-link groups.

      @param[in] pfi The PeptidoformIon to serialize
      @return JSON string representation of the PeptidoformIon
    */
    static String toJSON(const PeptidoformIon& pfi);

    /**
      @brief Construct a PeptidoformIon from JSON string

      Deserializes a JSON string back into a PeptidoformIon AST.

      @param[in] json_str JSON string representation of a PeptidoformIon
      @return The deserialized PeptidoformIon
      @throws Exception::ParseError if the JSON is malformed or missing required fields
    */
    static PeptidoformIon peptidoformIonFromJSON(const String& json_str);


    //==========================================================================
    // Static Methods - AASequence Conversion
    //==========================================================================

    /**
      @brief Resolve all modifications in a Peptidoform using ModificationsDB

      Looks up each modification tag (CV accession, named mod, mass delta) in
      ModificationsDB and stores the resolved ResidueModification pointer.

      @param[in,out] pf The Peptidoform to resolve (modified in place)
      @note Modifications that cannot be resolved will have resolved_mod = nullptr
    */
    static void resolveModifications(Peptidoform& pf);

    /**
      @brief Convert a Peptidoform to an OpenMS AASequence

      @param[in] pf The Peptidoform to convert
      @param[in] policy How to handle unconvertible modifications
      @return The equivalent AASequence
      @throws Exception::ConversionError if STRICT policy and conversion not possible

      @note Call resolveModifications() first, or this will resolve automatically
    */
    static AASequence toAASequence(
      const Peptidoform& pf,
      ConversionPolicy policy = ConversionPolicy::FAIL_ON_LOSS);

    /**
      @brief Create a Peptidoform from an OpenMS AASequence

      Converts an AASequence with modifications to ProForma notation.
      Uses CV accessions (UNIMOD) where available, otherwise named modifications.

      @param[in] seq The AASequence to convert
      @return The equivalent Peptidoform AST
    */
    static Peptidoform fromAASequence(const AASequence& seq);

    /**
      @brief Check if a Peptidoform can be fully represented as an AASequence

      Returns true if all modifications can be resolved and there are no
      unsupported features (ambiguous regions, cross-links, etc.)

      @param[in] pf The Peptidoform to check
      @return True if conversion is possible without issues
    */
    static bool isRepresentableAsAASequence(const Peptidoform& pf);

    /**
      @brief Get a list of all issues that would arise during AASequence conversion

      Returns detailed information about every aspect of the Peptidoform that
      cannot be represented in an AASequence.

      @param[in] pf The Peptidoform to analyze
      @return Vector of conversion issues (empty if fully convertible)
    */
    static std::vector<ConversionIssue> getAASequenceConversionIssues(const Peptidoform& pf);


    //==========================================================================
    // Static Methods - Mass Calculation
    //==========================================================================

    /**
      @brief Check if mass can be calculated for a Peptidoform

      Returns true if all components have known masses:
      - All amino acids are standard residues
      - All modifications are either resolved or have explicit mass deltas
      - No ambiguous regions with different possible amino acids

      @param[in] pf The Peptidoform to check (modifications will be resolved if needed)
      @return True if mass calculation is possible
    */
    static bool canCalculateMass(const Peptidoform& pf);

    /**
      @brief Check if mass can be calculated for a PeptidoformIon

      Returns true if mass can be calculated for all chains.
      Cross-links are handled correctly (cross-linker mass counted once).

      @param[in] pfi The PeptidoformIon to check
      @return True if mass calculation is possible
    */
    static bool canCalculateMass(const PeptidoformIon& pfi);

    /**
      @brief Get issues preventing mass calculation for a Peptidoform

      Returns detailed information about components that prevent mass calculation.

      @param[in] pf The Peptidoform to analyze
      @return Vector of issues (empty if mass can be calculated)
    */
    static std::vector<ConversionIssue> getMassCalculationIssues(const Peptidoform& pf);

    /**
      @brief Get issues preventing mass calculation for a PeptidoformIon

      Returns detailed information about components that prevent mass calculation
      across all chains.

      @param[in] pfi The PeptidoformIon to analyze
      @return Vector of issues (empty if mass can be calculated)
    */
    static std::vector<ConversionIssue> getMassCalculationIssues(const PeptidoformIon& pfi);

    /**
      @brief Calculate monoisotopic mass of a Peptidoform

      Calculates the neutral monoisotopic mass including:
      - Amino acid residue masses
      - Terminal H2O mass
      - All modification mass deltas
      - Unlocalised and labile modifications (included in total)
      - Global modifications applied to matching residues

      @param[in] pf The Peptidoform to calculate mass for
      @return Monoisotopic mass in Daltons
      @throws Exception::InvalidValue if mass cannot be calculated (use canCalculateMass() first)

      @note Modifications are resolved automatically if not already resolved
    */
    static double getMonoWeight(const Peptidoform& pf);

    /**
      @brief Calculate monoisotopic mass of a PeptidoformIon

      For cross-linked peptides, calculates the combined mass of all chains.
      Cross-linker masses are counted only once per cross-link group.

      @param[in] pfi The PeptidoformIon to calculate mass for
      @return Monoisotopic mass in Daltons
      @throws Exception::InvalidValue if mass cannot be calculated
      @throws Exception::InvalidValue if pfi is chimeric (use getMonoWeight on individual chains)
    */
    static double getMonoWeight(const PeptidoformIon& pfi);

    /**
      @brief Calculate m/z for a PeptidoformIon at its specified charge state

      Uses the charge state from the PeptidoformIon if present.

      @param[in] pfi The PeptidoformIon with charge state
      @return m/z value
      @throws Exception::InvalidValue if mass cannot be calculated or no charge state
    */
    static double getMZ(const PeptidoformIon& pfi);

    /**
      @brief Calculate m/z for a Peptidoform at a given charge state

      @param[in] pf The Peptidoform to calculate m/z for
      @param[in] charge The charge state (must be non-zero)
      @return m/z value
      @throws Exception::InvalidValue if mass cannot be calculated or charge is zero
    */
    static double getMZ(const Peptidoform& pf, int charge);

    // ---- Non-throwing variants (single-pass, efficient) ----

    /**
      @brief Try to calculate monoisotopic mass of a Peptidoform (non-throwing)

      Single-pass calculation that resolves modifications and calculates mass.
      More efficient than calling canCalculateMass() followed by getMonoWeight().

      @param[in] pf The Peptidoform to calculate mass for
      @return Monoisotopic mass in Daltons, or std::nullopt if calculation not possible
    */
    static std::optional<double> tryGetMonoWeight(const Peptidoform& pf);

    /**
      @brief Try to calculate monoisotopic mass with diagnostic information

      Single-pass calculation that also collects any issues preventing calculation.

      @param[in] pf The Peptidoform to calculate mass for
      @param[out] issues_out Vector to receive any issues (cleared first)
      @return Monoisotopic mass in Daltons, or std::nullopt if calculation not possible
    */
    static std::optional<double> tryGetMonoWeight(const Peptidoform& pf,
                                                   std::vector<ConversionIssue>& issues_out);

    /**
      @brief Try to calculate monoisotopic mass of a PeptidoformIon (non-throwing)

      @param[in] pfi The PeptidoformIon to calculate mass for
      @return Monoisotopic mass in Daltons, or std::nullopt if calculation not possible
      @note Returns std::nullopt for chimeric spectra (calculate per-chain instead)
    */
    static std::optional<double> tryGetMonoWeight(const PeptidoformIon& pfi);

    /**
      @brief Try to calculate monoisotopic mass of PeptidoformIon with diagnostics

      @param[in] pfi The PeptidoformIon to calculate mass for
      @param[out] issues_out Vector to receive any issues (cleared first)
      @return Monoisotopic mass in Daltons, or std::nullopt if calculation not possible
    */
    static std::optional<double> tryGetMonoWeight(const PeptidoformIon& pfi,
                                                   std::vector<ConversionIssue>& issues_out);

    /**
      @brief Try to calculate m/z for a Peptidoform (non-throwing)

      @param[in] pf The Peptidoform to calculate m/z for
      @param[in] charge The charge state (must be non-zero)
      @return m/z value, or std::nullopt if calculation not possible or charge is zero
    */
    static std::optional<double> tryGetMZ(const Peptidoform& pf, int charge);

    /**
      @brief Try to calculate m/z for a Peptidoform with diagnostics

      @param[in] pf The Peptidoform to calculate m/z for
      @param[in] charge The charge state (must be non-zero)
      @param[out] issues_out Vector to receive any issues (cleared first)
      @return m/z value, or std::nullopt if calculation not possible or charge is zero
    */
    static std::optional<double> tryGetMZ(const Peptidoform& pf, int charge,
                                           std::vector<ConversionIssue>& issues_out);

    /**
      @brief Try to calculate m/z for a PeptidoformIon (non-throwing)

      @param[in] pfi The PeptidoformIon with charge state
      @return m/z value, or std::nullopt if calculation not possible or no charge
    */
    static std::optional<double> tryGetMZ(const PeptidoformIon& pfi);

    /**
      @brief Try to calculate m/z for a PeptidoformIon with diagnostics

      @param[in] pfi The PeptidoformIon with charge state
      @param[out] issues_out Vector to receive any issues (cleared first)
      @return m/z value, or std::nullopt if calculation not possible or no charge
    */
    static std::optional<double> tryGetMZ(const PeptidoformIon& pfi,
                                           std::vector<ConversionIssue>& issues_out);


    //==========================================================================
    // Static Methods - Spectrum Generation
    //==========================================================================

    /**
      @brief Check if a theoretical spectrum can be generated for a Peptidoform

      Returns true if the Peptidoform can be converted to AASequence and
      fragmented. Use getSpectrumGenerationIssues() to get detailed diagnostics.

      @param[in] pf The Peptidoform to check
      @return True if spectrum generation is possible
    */
    static bool canGenerateSpectrum(const Peptidoform& pf);

    /**
      @brief Check if a theoretical spectrum can be generated for a PeptidoformIon

      Returns true if the PeptidoformIon can be fragmented. For cross-linked
      peptides, both chains must be convertible. Chimeric spectra are not supported.

      @param[in] pfi The PeptidoformIon to check
      @return True if spectrum generation is possible
    */
    static bool canGenerateSpectrum(const PeptidoformIon& pfi);

    /**
      @brief Get issues preventing spectrum generation for a Peptidoform

      @param[in] pf The Peptidoform to analyze
      @return Vector of issues (empty if spectrum can be generated)
    */
    static std::vector<ConversionIssue> getSpectrumGenerationIssues(const Peptidoform& pf);

    /**
      @brief Get issues preventing spectrum generation for a PeptidoformIon

      @param[in] pfi The PeptidoformIon to analyze
      @return Vector of issues (empty if spectrum can be generated)
    */
    static std::vector<ConversionIssue> getSpectrumGenerationIssues(const PeptidoformIon& pfi);

    /**
      @brief Generate a theoretical MS/MS spectrum for a Peptidoform

      Converts the Peptidoform to AASequence and uses TheoreticalSpectrumGenerator
      to generate fragment ions based on the specified ion types.

      @param[in] pf The Peptidoform to fragment
      @param[in] min_charge Minimum fragment ion charge state
      @param[in] max_charge Maximum fragment ion charge state
      @param[in] ion_types String specifying which ion types to generate:
                           'a','b','c','x','y','z' for ion series,
                           'M' for precursor peaks, 'I' for immonium ions.
                           Example: "by" for b/y ions, "byM" for b/y + precursor
      @param[in] add_losses If true, include neutral loss peaks (H2O, NH3)
      @param[in] add_metainfo If true, include ion annotations in spectrum
      @return MSSpectrum with theoretical peaks
      @throws Exception::ConversionError if spectrum generation fails
    */
    static MSSpectrum generateSpectrum(
      const Peptidoform& pf,
      int min_charge = 1,
      int max_charge = 1,
      const std::string& ion_types = "by",
      bool add_losses = false,
      bool add_metainfo = true);

    /**
      @brief Generate a theoretical MS/MS spectrum for a PeptidoformIon

      For single-chain peptides, uses TheoreticalSpectrumGenerator.
      For cross-linked peptides (// separator), uses TheoreticalSpectrumGeneratorXLMS.
      Chimeric spectra are not supported.

      @param[in] pfi The PeptidoformIon to fragment
      @param[in] min_charge Minimum fragment ion charge state
      @param[in] max_charge Maximum fragment ion charge state
      @param[in] ion_types String specifying which ion types to generate:
                           'a','b','c','x','y','z' for ion series,
                           'M' for precursor peaks, 'I' for immonium ions.
                           Example: "by" for b/y ions, "abyM" for a/b/y + precursor
      @param[in] add_losses If true, include neutral loss peaks
      @param[in] add_metainfo If true, include ion annotations
      @return MSSpectrum with theoretical peaks
      @throws Exception::ConversionError if spectrum generation fails
    */
    static MSSpectrum generateSpectrum(
      const PeptidoformIon& pfi,
      int min_charge = 1,
      int max_charge = 1,
      const std::string& ion_types = "by",
      bool add_losses = false,
      bool add_metainfo = true);


    //==========================================================================
    // Utility Functions
    //==========================================================================

    /**
      @brief Convert error code to human-readable string.

      @param[in] code The error code to convert.
      @return A human-readable string describing the error code.
    */
    static const char* errorCodeToString(ErrorCode code);
  }; // class ProForma

} // namespace OpenMS
