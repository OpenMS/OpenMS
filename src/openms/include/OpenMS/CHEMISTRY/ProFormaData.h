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

#include <nlohmann/json.hpp>

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
    STRICT,            ///< Fail if any modification cannot be fully represented
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
    - Connect cross-linked sites: [XLMOD:02001#XL1]...[#XL1]
    - Mark branch points: [#BRANCH]
    - Group ambiguously localized modifications: [Phospho#g1(0.90)]

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI Label
  {
    /// The type of label
    enum class Type
    {
      CROSSLINK,  ///< Cross-link label (e.g., #XL1)
      BRANCH,     ///< Branch point label (#BRANCH)
      AMBIGUOUS   ///< Ambiguous localization group (e.g., #g1)
    };

    Type type;                        ///< The label type
    String identifier;                ///< The label identifier (e.g., "XL1", "BRANCH", "g1")
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

    ProForma notation: <[TMT6plex]@K,N-term>

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

    ProForma notation: <13C> or <15N> or <D>

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI IsotopeReplacement
  {
    String isotope;  ///< The isotope specification (e.g., "13C", "15N", "D")
  };


  /**
    @brief Variant type for global modification entries

    A GlobalModEntry can be either:
    - IsotopeReplacement: Global isotope substitution (<13C>)
    - GlobalModification: Position-specific global mod (<[TMT6plex]@K>)

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
    - Global modifications (<13C>, <[TMT6plex]@K>)
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
    std::vector<GlobalModEntry> global_mods;          ///< Global modifications: <13C>, <[TMT6plex]@K>
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

    Derived during parsing from matching #XL labels.

    @ingroup Chemistry
  */
  struct OPENMS_DLLAPI CrossLinkGroup
  {
    String label;                                     ///< The cross-link label (e.g., "XL1")
    std::vector<std::pair<size_t, size_t>> sites;     ///< (chain_index, site_index) pairs
  };


  //--------------------------------------------------------------------------
  // JSON Serialization (nlohmann/json ADL pattern)
  //--------------------------------------------------------------------------

  /// @name CvDatabase enum JSON serialization
  /// @{

  /// Convert CvDatabase enum to JSON string
  inline void to_json(nlohmann::json& j, const CvDatabase& db)
  {
    switch (db)
    {
      case CvDatabase::UNIMOD: j = "UNIMOD"; break;
      case CvDatabase::MOD:    j = "MOD"; break;
      case CvDatabase::RESID:  j = "RESID"; break;
      case CvDatabase::XLMOD:  j = "XLMOD"; break;
      case CvDatabase::GNO:    j = "GNO"; break;
    }
  }

  /// Construct CvDatabase enum from JSON string
  inline void from_json(const nlohmann::json& j, CvDatabase& db)
  {
    std::string s = j.get<std::string>();
    if (s == "UNIMOD") db = CvDatabase::UNIMOD;
    else if (s == "MOD") db = CvDatabase::MOD;
    else if (s == "RESID") db = CvDatabase::RESID;
    else if (s == "XLMOD") db = CvDatabase::XLMOD;
    else if (s == "GNO") db = CvDatabase::GNO;
    else throw std::invalid_argument("Unknown CvDatabase: " + s);
  }

  /// @}


  /// @name CvAccession JSON serialization
  /// @{

  /// Convert CvAccession to JSON object
  inline void to_json(nlohmann::json& j, const CvAccession& cv)
  {
    j = nlohmann::json{{"database", cv.database}, {"accession", static_cast<std::string>(cv.accession)}};
  }

  /// Construct CvAccession from JSON object
  inline void from_json(const nlohmann::json& j, CvAccession& cv)
  {
    j.at("database").get_to(cv.database);
    cv.accession = j.at("accession").get<std::string>();
  }

  /// @}


  /// @name NamedMod JSON serialization
  /// @{

  /// Convert NamedMod to JSON object
  inline void to_json(nlohmann::json& j, const NamedMod& nm)
  {
    j = nlohmann::json{{"name", static_cast<std::string>(nm.name)}};
    if (nm.cv_hint.has_value())
    {
      j["cv_hint"] = nm.cv_hint.value();
    }
  }

  /// Construct NamedMod from JSON object
  inline void from_json(const nlohmann::json& j, NamedMod& nm)
  {
    nm.name = j.at("name").get<std::string>();
    if (j.contains("cv_hint") && !j.at("cv_hint").is_null())
    {
      nm.cv_hint = j.at("cv_hint").get<CvDatabase>();
    }
    else
    {
      nm.cv_hint = std::nullopt;
    }
  }

  /// @}


  /// @name MassDelta::Source enum JSON serialization
  /// @{

  /// Convert MassDelta::Source enum to JSON string
  inline void to_json(nlohmann::json& j, const MassDelta::Source& src)
  {
    switch (src)
    {
      case MassDelta::Source::NONE: j = "NONE"; break;
      case MassDelta::Source::OBS:  j = "OBS"; break;
      case MassDelta::Source::U:    j = "U"; break;
      case MassDelta::Source::M:    j = "M"; break;
      case MassDelta::Source::R:    j = "R"; break;
      case MassDelta::Source::X:    j = "X"; break;
      case MassDelta::Source::G:    j = "G"; break;
    }
  }

  /// Construct MassDelta::Source enum from JSON string
  inline void from_json(const nlohmann::json& j, MassDelta::Source& src)
  {
    std::string s = j.get<std::string>();
    if (s == "NONE") src = MassDelta::Source::NONE;
    else if (s == "OBS") src = MassDelta::Source::OBS;
    else if (s == "U") src = MassDelta::Source::U;
    else if (s == "M") src = MassDelta::Source::M;
    else if (s == "R") src = MassDelta::Source::R;
    else if (s == "X") src = MassDelta::Source::X;
    else if (s == "G") src = MassDelta::Source::G;
    else throw std::invalid_argument("Unknown MassDelta::Source: " + s);
  }

  /// @}


  /// @name MassDelta JSON serialization
  /// @{

  /// Convert MassDelta to JSON object
  inline void to_json(nlohmann::json& j, const MassDelta& md)
  {
    j = nlohmann::json{
      {"source", md.source},
      {"mass", md.mass},
      {"original_text", static_cast<std::string>(md.original_text)}
    };
  }

  /// Construct MassDelta from JSON object
  inline void from_json(const nlohmann::json& j, MassDelta& md)
  {
    j.at("source").get_to(md.source);
    j.at("mass").get_to(md.mass);
    md.original_text = j.at("original_text").get<std::string>();
  }

  /// @}


  /// @name FormulaTag JSON serialization
  /// @{

  /// Convert FormulaTag to JSON object
  inline void to_json(nlohmann::json& j, const FormulaTag& ft)
  {
    j = nlohmann::json{{"formula_string", static_cast<std::string>(ft.formula_string)}};
    if (ft.charge.has_value())
    {
      j["charge"] = ft.charge.value();
    }
  }

  /// Construct FormulaTag from JSON object
  inline void from_json(const nlohmann::json& j, FormulaTag& ft)
  {
    ft.formula_string = j.at("formula_string").get<std::string>();
    if (j.contains("charge") && !j.at("charge").is_null())
    {
      ft.charge = j.at("charge").get<int>();
    }
    else
    {
      ft.charge = std::nullopt;
    }
  }

  /// @}


  /// @name GlycanComposition::Monosaccharide JSON serialization
  /// @{

  /// Convert Monosaccharide variant to JSON object
  inline void to_json(nlohmann::json& j, const GlycanComposition::Monosaccharide& mono)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, String>)
      {
        j = nlohmann::json{{"type", "name"}, {"value", static_cast<std::string>(arg)}};
      }
      else if constexpr (std::is_same_v<T, FormulaTag>)
      {
        j = nlohmann::json{{"type", "formula"}, {"value", arg}};
      }
    }, mono);
  }

  /// Construct Monosaccharide variant from JSON object
  inline void from_json(const nlohmann::json& j, GlycanComposition::Monosaccharide& mono)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "name")
    {
      mono = String(j.at("value").get<std::string>());
    }
    else if (type == "formula")
    {
      mono = j.at("value").get<FormulaTag>();
    }
    else
    {
      throw std::invalid_argument("Unknown Monosaccharide type: " + type);
    }
  }

  /// @}


  /// @name GlycanComposition JSON serialization
  /// @{

  /// Convert GlycanComposition to JSON object
  inline void to_json(nlohmann::json& j, const GlycanComposition& gc)
  {
    j = nlohmann::json::array();
    for (const auto& [mono, count] : gc.components)
    {
      nlohmann::json component;
      component["monosaccharide"] = mono;
      component["count"] = count;
      j.push_back(component);
    }
  }

  /// Construct GlycanComposition from JSON object
  inline void from_json(const nlohmann::json& j, GlycanComposition& gc)
  {
    gc.components.clear();
    for (const auto& item : j)
    {
      GlycanComposition::Monosaccharide mono;
      from_json(item.at("monosaccharide"), mono);
      int count = item.at("count").get<int>();
      gc.components.emplace_back(mono, count);
    }
  }

  /// @}


  /// @name InfoTag JSON serialization
  /// @{

  /// Convert InfoTag to JSON object
  inline void to_json(nlohmann::json& j, const InfoTag& it)
  {
    j = nlohmann::json{{"text", static_cast<std::string>(it.text)}};
  }

  /// Construct InfoTag from JSON object
  inline void from_json(const nlohmann::json& j, InfoTag& it)
  {
    it.text = j.at("text").get<std::string>();
  }

  /// @}


  /// @name PositionConstraint JSON serialization
  /// @{

  /// Convert PositionConstraint to JSON object
  inline void to_json(nlohmann::json& j, const PositionConstraint& pc)
  {
    std::string residue_str;
    for (char c : pc.residues)
    {
      residue_str += c;
    }
    j = nlohmann::json{{"residues", residue_str}};
  }

  /// Construct PositionConstraint from JSON object
  inline void from_json(const nlohmann::json& j, PositionConstraint& pc)
  {
    std::string residue_str = j.at("residues").get<std::string>();
    pc.residues.clear();
    for (char c : residue_str)
    {
      pc.residues.push_back(c);
    }
  }

  /// @}


  /// @name ModificationTag variant JSON serialization
  /// @{

  /// Convert ModificationTag variant to JSON object
  inline void to_json(nlohmann::json& j, const ModificationTag& tag)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, CvAccession>)
      {
        j = nlohmann::json{{"type", "cv_accession"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, NamedMod>)
      {
        j = nlohmann::json{{"type", "named_mod"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, MassDelta>)
      {
        j = nlohmann::json{{"type", "mass_delta"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, FormulaTag>)
      {
        j = nlohmann::json{{"type", "formula"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, GlycanComposition>)
      {
        j = nlohmann::json{{"type", "glycan"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, InfoTag>)
      {
        j = nlohmann::json{{"type", "info"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, PositionConstraint>)
      {
        j = nlohmann::json{{"type", "position"}, {"value", arg}};
      }
    }, tag);
  }

  /// Construct ModificationTag variant from JSON object
  inline void from_json(const nlohmann::json& j, ModificationTag& tag)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "cv_accession")
    {
      tag = j.at("value").get<CvAccession>();
    }
    else if (type == "named_mod")
    {
      tag = j.at("value").get<NamedMod>();
    }
    else if (type == "mass_delta")
    {
      tag = j.at("value").get<MassDelta>();
    }
    else if (type == "formula")
    {
      tag = j.at("value").get<FormulaTag>();
    }
    else if (type == "glycan")
    {
      tag = j.at("value").get<GlycanComposition>();
    }
    else if (type == "info")
    {
      tag = j.at("value").get<InfoTag>();
    }
    else if (type == "position")
    {
      tag = j.at("value").get<PositionConstraint>();
    }
    else
    {
      throw std::invalid_argument("Unknown ModificationTag type: " + type);
    }
  }

  /// @}


  /// @name Label::Type enum JSON serialization
  /// @{

  /// Convert Label::Type enum to JSON string
  inline void to_json(nlohmann::json& j, const Label::Type& lt)
  {
    switch (lt)
    {
      case Label::Type::CROSSLINK: j = "CROSSLINK"; break;
      case Label::Type::BRANCH:    j = "BRANCH"; break;
      case Label::Type::AMBIGUOUS: j = "AMBIGUOUS"; break;
    }
  }

  /// Construct Label::Type enum from JSON string
  inline void from_json(const nlohmann::json& j, Label::Type& lt)
  {
    std::string s = j.get<std::string>();
    if (s == "CROSSLINK") lt = Label::Type::CROSSLINK;
    else if (s == "BRANCH") lt = Label::Type::BRANCH;
    else if (s == "AMBIGUOUS") lt = Label::Type::AMBIGUOUS;
    else throw std::invalid_argument("Unknown Label::Type: " + s);
  }

  /// @}


  /// @name Label JSON serialization
  /// @{

  /// Convert Label to JSON object
  inline void to_json(nlohmann::json& j, const Label& lbl)
  {
    j = nlohmann::json{
      {"type", lbl.type},
      {"identifier", static_cast<std::string>(lbl.identifier)}
    };
    if (lbl.score.has_value())
    {
      j["score"] = lbl.score.value();
    }
  }

  /// Construct Label from JSON object
  inline void from_json(const nlohmann::json& j, Label& lbl)
  {
    j.at("type").get_to(lbl.type);
    lbl.identifier = j.at("identifier").get<std::string>();
    if (j.contains("score") && !j.at("score").is_null())
    {
      lbl.score = j.at("score").get<double>();
    }
    else
    {
      lbl.score = std::nullopt;
    }
  }

  /// @}


  /// @name Modification JSON serialization
  /// @{

  /// Convert Modification to JSON object
  inline void to_json(nlohmann::json& j, const Modification& mod)
  {
    j = nlohmann::json::array();
    for (const auto& [tag, label] : mod.alternatives)
    {
      nlohmann::json alt;
      alt["tag"] = tag;
      if (label.has_value())
      {
        alt["label"] = label.value();
      }
      j.push_back(alt);
    }
  }

  /// Construct Modification from JSON object
  inline void from_json(const nlohmann::json& j, Modification& mod)
  {
    mod.alternatives.clear();
    for (const auto& item : j)
    {
      ModificationTag tag = item.at("tag").get<ModificationTag>();
      std::optional<Label> label;
      if (item.contains("label") && !item.at("label").is_null())
      {
        label = item.at("label").get<Label>();
      }
      mod.alternatives.emplace_back(tag, label);
    }
  }

  /// @}


  /// @name SequenceElement JSON serialization
  /// @{

  /// Convert SequenceElement to JSON object
  inline void to_json(nlohmann::json& j, const SequenceElement& se)
  {
    j = nlohmann::json{
      {"amino_acid", std::string(1, se.amino_acid)},
      {"modifications", se.modifications}
    };
  }

  /// Construct SequenceElement from JSON object
  inline void from_json(const nlohmann::json& j, SequenceElement& se)
  {
    std::string aa = j.at("amino_acid").get<std::string>();
    if (aa.empty())
    {
      throw std::invalid_argument("amino_acid cannot be empty");
    }
    se.amino_acid = aa[0];
    j.at("modifications").get_to(se.modifications);
  }

  /// @}


  /// @name AmbiguousRegion JSON serialization
  /// @{

  /// Convert AmbiguousRegion to JSON object
  inline void to_json(nlohmann::json& j, const AmbiguousRegion& ar)
  {
    j = nlohmann::json{{"elements", ar.elements}};
  }

  /// Construct AmbiguousRegion from JSON object
  inline void from_json(const nlohmann::json& j, AmbiguousRegion& ar)
  {
    j.at("elements").get_to(ar.elements);
  }

  /// @}


  /// @name ModifiedRange JSON serialization
  /// @{

  /// Convert ModifiedRange to JSON object
  inline void to_json(nlohmann::json& j, const ModifiedRange& mr)
  {
    j = nlohmann::json{
      {"elements", mr.elements},
      {"modifications", mr.modifications}
    };
  }

  /// Construct ModifiedRange from JSON object
  inline void from_json(const nlohmann::json& j, ModifiedRange& mr)
  {
    j.at("elements").get_to(mr.elements);
    j.at("modifications").get_to(mr.modifications);
  }

  /// @}


  /// @name SequenceSection variant JSON serialization
  /// @{

  /// Convert SequenceSection variant to JSON object
  inline void to_json(nlohmann::json& j, const SequenceSection& ss)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, SequenceElement>)
      {
        j = nlohmann::json{{"type", "element"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, AmbiguousRegion>)
      {
        j = nlohmann::json{{"type", "ambiguous_region"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, ModifiedRange>)
      {
        j = nlohmann::json{{"type", "modified_range"}, {"value", arg}};
      }
    }, ss);
  }

  /// Construct SequenceSection variant from JSON object
  inline void from_json(const nlohmann::json& j, SequenceSection& ss)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "element")
    {
      ss = j.at("value").get<SequenceElement>();
    }
    else if (type == "ambiguous_region")
    {
      ss = j.at("value").get<AmbiguousRegion>();
    }
    else if (type == "modified_range")
    {
      ss = j.at("value").get<ModifiedRange>();
    }
    else
    {
      throw std::invalid_argument("Unknown SequenceSection type: " + type);
    }
  }

  /// @}


  /// @name UnlocalisedMod JSON serialization
  /// @{

  /// Convert UnlocalisedMod to JSON object
  inline void to_json(nlohmann::json& j, const UnlocalisedMod& um)
  {
    j = nlohmann::json{{"modifications", um.modifications}};
    if (um.occurrence.has_value())
    {
      j["occurrence"] = um.occurrence.value();
    }
  }

  /// Construct UnlocalisedMod from JSON object
  inline void from_json(const nlohmann::json& j, UnlocalisedMod& um)
  {
    j.at("modifications").get_to(um.modifications);
    if (j.contains("occurrence") && !j.at("occurrence").is_null())
    {
      um.occurrence = j.at("occurrence").get<int>();
    }
    else
    {
      um.occurrence = std::nullopt;
    }
  }

  /// @}


  /// @name LabileModification JSON serialization
  /// @{

  /// Convert LabileModification to JSON object
  inline void to_json(nlohmann::json& j, const LabileModification& lm)
  {
    j = nlohmann::json{{"modification", lm.modification}};
  }

  /// Construct LabileModification from JSON object
  inline void from_json(const nlohmann::json& j, LabileModification& lm)
  {
    j.at("modification").get_to(lm.modification);
  }

  /// @}


  /// @name GlobalModification JSON serialization
  /// @{

  /// Convert GlobalModification to JSON object
  inline void to_json(nlohmann::json& j, const GlobalModification& gm)
  {
    std::vector<std::string> locs;
    for (const auto& loc : gm.locations)
    {
      locs.push_back(static_cast<std::string>(loc));
    }
    j = nlohmann::json{
      {"modification", gm.modification},
      {"locations", locs}
    };
  }

  /// Construct GlobalModification from JSON object
  inline void from_json(const nlohmann::json& j, GlobalModification& gm)
  {
    j.at("modification").get_to(gm.modification);
    gm.locations.clear();
    for (const auto& loc : j.at("locations"))
    {
      gm.locations.push_back(String(loc.get<std::string>()));
    }
  }

  /// @}


  /// @name IsotopeReplacement JSON serialization
  /// @{

  /// Convert IsotopeReplacement to JSON object
  inline void to_json(nlohmann::json& j, const IsotopeReplacement& ir)
  {
    j = nlohmann::json{{"isotope", static_cast<std::string>(ir.isotope)}};
  }

  /// Construct IsotopeReplacement from JSON object
  inline void from_json(const nlohmann::json& j, IsotopeReplacement& ir)
  {
    ir.isotope = j.at("isotope").get<std::string>();
  }

  /// @}


  /// @name GlobalModEntry variant JSON serialization
  /// @{

  /// Convert GlobalModEntry variant to JSON object
  inline void to_json(nlohmann::json& j, const GlobalModEntry& gme)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, IsotopeReplacement>)
      {
        j = nlohmann::json{{"type", "isotope_replacement"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, GlobalModification>)
      {
        j = nlohmann::json{{"type", "global_modification"}, {"value", arg}};
      }
    }, gme);
  }

  /// Construct GlobalModEntry variant from JSON object
  inline void from_json(const nlohmann::json& j, GlobalModEntry& gme)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "isotope_replacement")
    {
      gme = j.at("value").get<IsotopeReplacement>();
    }
    else if (type == "global_modification")
    {
      gme = j.at("value").get<GlobalModification>();
    }
    else
    {
      throw std::invalid_argument("Unknown GlobalModEntry type: " + type);
    }
  }

  /// @}


  /// @name AdductIon JSON serialization
  /// @{

  /// Convert AdductIon to JSON object
  inline void to_json(nlohmann::json& j, const AdductIon& ai)
  {
    j = nlohmann::json{
      {"formula", static_cast<std::string>(ai.formula)},
      {"charge", ai.charge}
    };
    if (ai.occurrence.has_value())
    {
      j["occurrence"] = ai.occurrence.value();
    }
  }

  /// Construct AdductIon from JSON object
  inline void from_json(const nlohmann::json& j, AdductIon& ai)
  {
    ai.formula = j.at("formula").get<std::string>();
    j.at("charge").get_to(ai.charge);
    if (j.contains("occurrence") && !j.at("occurrence").is_null())
    {
      ai.occurrence = j.at("occurrence").get<int>();
    }
    else
    {
      ai.occurrence = std::nullopt;
    }
  }

  /// @}


  /// @name ChargeState variant JSON serialization
  /// @{

  /// Convert ChargeState variant to JSON object
  inline void to_json(nlohmann::json& j, const ChargeState& cs)
  {
    std::visit([&j](auto&& arg) {
      using T = std::decay_t<decltype(arg)>;
      if constexpr (std::is_same_v<T, int>)
      {
        j = nlohmann::json{{"type", "simple"}, {"value", arg}};
      }
      else if constexpr (std::is_same_v<T, std::vector<AdductIon>>)
      {
        j = nlohmann::json{{"type", "adducts"}, {"value", arg}};
      }
    }, cs);
  }

  /// Construct ChargeState variant from JSON object
  inline void from_json(const nlohmann::json& j, ChargeState& cs)
  {
    std::string type = j.at("type").get<std::string>();
    if (type == "simple")
    {
      cs = j.at("value").get<int>();
    }
    else if (type == "adducts")
    {
      cs = j.at("value").get<std::vector<AdductIon>>();
    }
    else
    {
      throw std::invalid_argument("Unknown ChargeState type: " + type);
    }
  }

  /// @}


  /// @name Peptidoform JSON serialization
  /// @{

  /// Convert Peptidoform to JSON object
  inline void to_json(nlohmann::json& j, const Peptidoform& pf)
  {
    j = nlohmann::json{
      {"global_mods", pf.global_mods},
      {"unlocalised_mods", pf.unlocalised_mods},
      {"labile_mods", pf.labile_mods},
      {"n_term_mods", pf.n_term_mods},
      {"sequence", pf.sequence},
      {"c_term_mods", pf.c_term_mods}
    };
    if (pf.name.has_value())
    {
      j["name"] = static_cast<std::string>(pf.name.value());
    }
    if (pf.charge.has_value())
    {
      j["charge"] = pf.charge.value();
    }
  }

  /// Construct Peptidoform from JSON object
  inline void from_json(const nlohmann::json& j, Peptidoform& pf)
  {
    if (j.contains("global_mods"))
    {
      j.at("global_mods").get_to(pf.global_mods);
    }
    j.at("unlocalised_mods").get_to(pf.unlocalised_mods);
    j.at("labile_mods").get_to(pf.labile_mods);
    j.at("n_term_mods").get_to(pf.n_term_mods);
    j.at("sequence").get_to(pf.sequence);
    j.at("c_term_mods").get_to(pf.c_term_mods);
    if (j.contains("name") && !j.at("name").is_null())
    {
      pf.name = String(j.at("name").get<std::string>());
    }
    else
    {
      pf.name = std::nullopt;
    }
    if (j.contains("charge") && !j.at("charge").is_null())
    {
      pf.charge = j.at("charge").get<ChargeState>();
    }
    else
    {
      pf.charge = std::nullopt;
    }
  }

  /// @}


  /// @name PeptidoformIon JSON serialization
  /// @{

  /// Convert PeptidoformIon to JSON object
  inline void to_json(nlohmann::json& j, const PeptidoformIon& pfi)
  {
    j = nlohmann::json{{"chains", pfi.chains}, {"is_chimeric", pfi.is_chimeric}};
    if (pfi.name.has_value())
    {
      j["name"] = static_cast<std::string>(pfi.name.value());
    }
    if (pfi.charge.has_value())
    {
      j["charge"] = pfi.charge.value();
    }
  }

  /// Construct PeptidoformIon from JSON object
  inline void from_json(const nlohmann::json& j, PeptidoformIon& pfi)
  {
    j.at("chains").get_to(pfi.chains);
    if (j.contains("name") && !j.at("name").is_null())
    {
      pfi.name = String(j.at("name").get<std::string>());
    }
    else
    {
      pfi.name = std::nullopt;
    }
    if (j.contains("charge") && !j.at("charge").is_null())
    {
      pfi.charge = j.at("charge").get<ChargeState>();
    }
    else
    {
      pfi.charge = std::nullopt;
    }
    if (j.contains("is_chimeric"))
    {
      pfi.is_chimeric = j.at("is_chimeric").get<bool>();
    }
    else
    {
      pfi.is_chimeric = false;
    }
  }

  /// @}


  /// @name CrossLinkGroup JSON serialization
  /// @{

  /// Convert CrossLinkGroup to JSON object
  inline void to_json(nlohmann::json& j, const CrossLinkGroup& clg)
  {
    nlohmann::json sites_json = nlohmann::json::array();
    for (const auto& [chain_idx, site_idx] : clg.sites)
    {
      sites_json.push_back(nlohmann::json{{"chain_index", chain_idx}, {"site_index", site_idx}});
    }
    j = nlohmann::json{
      {"label", static_cast<std::string>(clg.label)},
      {"sites", sites_json}
    };
  }

  /// Construct CrossLinkGroup from JSON object
  inline void from_json(const nlohmann::json& j, CrossLinkGroup& clg)
  {
    clg.label = j.at("label").get<std::string>();
    clg.sites.clear();
    for (const auto& site : j.at("sites"))
    {
      size_t chain_idx = site.at("chain_index").get<size_t>();
      size_t site_idx = site.at("site_index").get<size_t>();
      clg.sites.emplace_back(chain_idx, site_idx);
    }
  }

  /// @}


  //--------------------------------------------------------------------------
  // Convenience methods for JSON string conversion
  //--------------------------------------------------------------------------

  /// Convert Peptidoform to JSON string
  inline String toJSON(const Peptidoform& pf)
  {
    nlohmann::json j = pf;
    return String(j.dump());
  }

  /// Construct Peptidoform from JSON string
  inline Peptidoform peptidoformFromJSON(const String& json_str)
  {
    nlohmann::json j = nlohmann::json::parse(static_cast<std::string>(json_str));
    return j.get<Peptidoform>();
  }

  /// Convert PeptidoformIon to JSON string
  inline String toJSON(const PeptidoformIon& pfi)
  {
    nlohmann::json j = pfi;
    return String(j.dump());
  }

  /// Construct PeptidoformIon from JSON string
  inline PeptidoformIon peptidoformIonFromJSON(const String& json_str)
  {
    nlohmann::json j = nlohmann::json::parse(static_cast<std::string>(json_str));
    return j.get<PeptidoformIon>();
  }

} // namespace OpenMS
