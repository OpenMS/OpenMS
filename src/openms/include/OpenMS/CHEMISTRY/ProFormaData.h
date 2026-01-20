// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>


#include <optional>
#include <utility>
#include <variant>
#include <vector>

namespace OpenMS
{

  // Forward declaration for resolved modification pointer
  class ResidueModification;

  /**
    @brief Conversion policy for transforming Peptidoform to AASequence

    Controls how the conversion handles modifications that cannot be directly
    represented in AASequence (e.g., unlocalised, labile, or ambiguous modifications).

    @ingroup Chemistry
  */
  enum class AASequenceConversionPolicy
  {
    STRICT_MODE,       ///< Fail if any modification cannot be fully represented
    DROP_UNLOCALISED,  ///< Drop unlocalised, labile, and global modifications
    BEST_EFFORT        ///< Try to convert as much as possible, skip unsupported
  };


  /**
    @brief Issue type for AASequence conversion problems

    @ingroup Chemistry
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
    @brief Description of a conversion issue from Peptidoform to AASequence

    Records problems encountered when attempting to convert a ProForma
    Peptidoform to an OpenMS AASequence representation.

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI ConversionIssue
  {
    ConversionIssueType type;   ///< The type of issue
    String description;         ///< Human-readable description
    size_t position;            ///< Position in sequence (SIZE_MAX if not position-specific)
  };


  /**
    @brief Write mode for ProForma string serialization

    Controls whether the output preserves original formatting (LOSSLESS) or
    produces a normalized, deterministic output (CANONICAL).

    @ingroup Chemistry
  */
  enum class ProFormaWriteMode
  {
    LOSSLESS,   ///< Preserve original spelling/formatting where possible (e.g., mass delta text)
    CANONICAL   ///< Normalized output: uppercase CV prefixes, sorted mods, 4 decimal places for masses
  };


  /**
    @brief Controlled vocabulary database prefix for modification accessions

    Identifies the source database for a modification accession in ProForma notation.
    Examples: UNIMOD:35, MOD:00046, XLMOD:02001, GNO:G59626AS

    @ingroup Chemistry
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
    @brief Controlled vocabulary accession for a modification

    Represents a modification specified by a CV accession number, e.g., UNIMOD:35 for Oxidation.
    The accession string contains only the identifier portion (e.g., "35" for UNIMOD:35).

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
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

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI CrossLinkGroup
  {
    String label;                                     ///< The cross-link label (e.g., `XL1`)
    std::vector<std::pair<size_t, size_t>> sites;     ///< (chain_index, site_index) pairs
  };


  //--------------------------------------------------------------------------
  // JSON serialization convenience functions (implementations in ProFormaDataJson.cpp)
  //--------------------------------------------------------------------------

  /// @name JSON serialization
  /// @{

  /**
    @brief Convert a Peptidoform to JSON string representation

    Serializes the complete Peptidoform AST including all modifications,
    terminal modifications, global modifications, and labels.

    @param[in] pf The Peptidoform to serialize
    @return JSON string representation of the Peptidoform
  */
  OPENMS_DLLAPI String toJSON(const Peptidoform& pf);

  /**
    @brief Construct a Peptidoform from JSON string

    Deserializes a JSON string back into a Peptidoform AST.

    @param[in] json_str JSON string representation of a Peptidoform
    @return The deserialized Peptidoform
    @throws Exception::ParseError if the JSON is malformed or missing required fields
  */
  OPENMS_DLLAPI Peptidoform peptidoformFromJSON(const String& json_str);

  /**
    @brief Convert a PeptidoformIon to JSON string representation

    Serializes the complete PeptidoformIon AST including all chains,
    charge state, adduct ions, and cross-link groups.

    @param[in] pfi The PeptidoformIon to serialize
    @return JSON string representation of the PeptidoformIon
  */
  OPENMS_DLLAPI String toJSON(const PeptidoformIon& pfi);

  /**
    @brief Construct a PeptidoformIon from JSON string

    Deserializes a JSON string back into a PeptidoformIon AST.

    @param[in] json_str JSON string representation of a PeptidoformIon
    @return The deserialized PeptidoformIon
    @throws Exception::ParseError if the JSON is malformed or missing required fields
  */
  OPENMS_DLLAPI PeptidoformIon peptidoformIonFromJSON(const String& json_str);

  /// @}


} // namespace OpenMS
