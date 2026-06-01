// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: OpenMS Team $
// $Authors: OpenMS Team $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>

#include <fstream>
#include <limits>
#include <map>
#include <vector>

namespace OpenMS
{
  /**
    @brief Represents a PEFF modification annotation.

    PEFF modifications can be in PSI-MOD, UNIMOD, or custom format.
    Position is 1-based, with 0 indicating unknown position (represented as ? in PEFF).
  */
  struct OPENMS_DLLAPI PEFFModification
  {
    Size position{0};         ///< 1-based position, 0 = unknown position (?)
    String accession;         ///< "MOD:00046", "UNIMOD:35", or custom
    String name;              ///< Human-readable name
    String optional_tag;      ///< Optional tag (last component of annotation tuple)
    UInt annotation_id{std::numeric_limits<UInt>::max()};  ///< Optional annotation identifier (when HasAnnotationIdentifiers=true), max() = not set

    enum class Type { PSI_MOD, UNIMOD, GENERIC };
    Type type{Type::GENERIC};

    PEFFModification() = default;
    PEFFModification(Size pos, const String& acc, const String& n, const String& tag = "", UInt aid = std::numeric_limits<UInt>::max())
      : position(pos), accession(acc), name(n), optional_tag(tag), annotation_id(aid)
    {
      if (accession.hasPrefix("MOD:"))
      {
        type = Type::PSI_MOD;
      }
      else if (accession.hasPrefix("UNIMOD:"))
      {
        type = Type::UNIMOD;
      }
    }

    bool operator==(const PEFFModification& rhs) const
    {
      return position == rhs.position && accession == rhs.accession &&
             name == rhs.name && type == rhs.type &&
             optional_tag == rhs.optional_tag && annotation_id == rhs.annotation_id;
    }
  };

  /**
    @brief Represents a simple PEFF variant (single amino acid substitution).

    Parsed from \\VariantSimple annotations.
  */
  struct OPENMS_DLLAPI PEFFVariantSimple
  {
    Size position{0};         ///< 1-based position
    char variant_aa{'\0'};    ///< Variant amino acid
    String optional_tag;      ///< Optional tag (last component of annotation tuple)
    UInt annotation_id{std::numeric_limits<UInt>::max()};  ///< Optional annotation identifier, max() = not set

    PEFFVariantSimple() = default;
    PEFFVariantSimple(Size pos, char aa, const String& tag = "", UInt aid = std::numeric_limits<UInt>::max())
      : position(pos), variant_aa(aa), optional_tag(tag), annotation_id(aid) {}

    bool operator==(const PEFFVariantSimple& rhs) const
    {
      return position == rhs.position && variant_aa == rhs.variant_aa &&
             optional_tag == rhs.optional_tag && annotation_id == rhs.annotation_id;
    }
  };

  /**
    @brief Represents a complex PEFF variant (insertion, deletion, or substitution of multiple amino acids).

    Parsed from \\VariantComplex annotations.
  */
  struct OPENMS_DLLAPI PEFFVariantComplex
  {
    Size start_position{0};   ///< 1-based start position
    Size end_position{0};     ///< 1-based end position
    String replacement;       ///< Replacement sequence (empty = deletion)
    String optional_tag;      ///< Optional tag (last component of annotation tuple)
    UInt annotation_id{std::numeric_limits<UInt>::max()};  ///< Optional annotation identifier, max() = not set

    PEFFVariantComplex() = default;
    PEFFVariantComplex(Size start, Size end, const String& repl, const String& tag = "", UInt aid = std::numeric_limits<UInt>::max())
      : start_position(start), end_position(end), replacement(repl), optional_tag(tag), annotation_id(aid) {}

    bool operator==(const PEFFVariantComplex& rhs) const
    {
      return start_position == rhs.start_position && end_position == rhs.end_position &&
             replacement == rhs.replacement && optional_tag == rhs.optional_tag &&
             annotation_id == rhs.annotation_id;
    }
  };

  /**
    @brief Represents a PEFF processed region (signal peptide, transit peptide, etc.).

    Parsed from \\Processed annotations.
  */
  struct OPENMS_DLLAPI PEFFProcessedRegion
  {
    Size start_position{0};   ///< 1-based start position
    Size end_position{0};     ///< 1-based end position
    String accession;         ///< PEFF CV accession (e.g., "PEFF:0001021")
    String name;              ///< Optional name (e.g., "signal peptide")
    String optional_tag;      ///< Optional tag (last component of annotation tuple)
    UInt annotation_id{std::numeric_limits<UInt>::max()};  ///< Optional annotation identifier, max() = not set

    PEFFProcessedRegion() = default;
    PEFFProcessedRegion(Size start, Size end, const String& acc, const String& n = "", const String& tag = "", UInt aid = std::numeric_limits<UInt>::max())
      : start_position(start), end_position(end), accession(acc), name(n), optional_tag(tag), annotation_id(aid) {}

    bool operator==(const PEFFProcessedRegion& rhs) const
    {
      return start_position == rhs.start_position && end_position == rhs.end_position &&
             accession == rhs.accession && name == rhs.name && optional_tag == rhs.optional_tag &&
             annotation_id == rhs.annotation_id;
    }
  };

  /**
    @brief Represents a disulfide bond annotation in PEFF.

    Parsed from \\DisulfideBond annotations.
  */
  struct OPENMS_DLLAPI PEFFDisulfideBond
  {
    String id1;                  ///< First cysteine reference (AnnotationIdentifier of the cysteine residue)
    String id2;                  ///< Second cysteine reference (AnnotationIdentifier of the cysteine residue)
    String optional_tag;         ///< Optional tag (e.g., "between chains")
    UInt annotation_id{std::numeric_limits<UInt>::max()};  ///< Optional annotation identifier, max() = not set

    PEFFDisulfideBond() = default;
    PEFFDisulfideBond(const String& i1, const String& i2, const String& tag = "", UInt aid = std::numeric_limits<UInt>::max())
      : id1(i1), id2(i2), optional_tag(tag), annotation_id(aid) {}

    bool operator==(const PEFFDisulfideBond& rhs) const
    {
      return id1 == rhs.id1 && id2 == rhs.id2 && optional_tag == rhs.optional_tag &&
             annotation_id == rhs.annotation_id;
    }
  };

  /**
    @brief Represents a single entry in a PEFF file with all annotations.

    Each entry corresponds to one description line and sequence in the PEFF file.
    The description line format per the PEFF spec is:

      >Prefix:DbUniqueId \\key=value \\key=value ...

    Where Prefix is the database prefix defined in the header block and DbUniqueId is
    the unique identifier within that database. The identifier field stores the full
    "Prefix:DbUniqueId" string, and the prefix field stores just the prefix portion.
  */
  struct OPENMS_DLLAPI PEFFEntry
  {
    // Basic fields
    String prefix;             ///< Database prefix from description line (e.g., "sp" from ">sp:P12345")
    String identifier;
    String sequence;

    // Metadata
    std::vector<String> protein_names;   ///< \\PName - may have multiple names
    String gene_name;                     ///< \\GName
    Int ncbi_tax_id{0};                  ///< \\NcbiTaxId or \\OX
    String taxonomy_name;                 ///< \\TaxName
    Size sequence_length{0};             ///< \\Length
    String sequence_version;              ///< \\SV
    String entry_version;                 ///< \\EV
    Int protein_existence{0};            ///< \\PE (1-5)
    String db_unique_id;                  ///< \\DbUniqueId
    String entry_id;                      ///< \\ID (e.g., NPM_HUMAN)
    std::vector<String> alt_accessions;   ///< \\AltAC - alternative accessions

    // Annotations
    std::vector<PEFFModification> modifications;
    std::vector<PEFFVariantSimple> simple_variants;
    std::vector<PEFFVariantComplex> complex_variants;
    std::vector<PEFFProcessedRegion> processed_regions;
    std::vector<PEFFDisulfideBond> disulfide_bonds;  ///< \\DisulfideBond
    std::vector<String> proteoforms;     ///< ProForma notation

    // Custom annotations
    std::map<String, String> custom_annotations;

    PEFFEntry() = default;

    PEFFEntry(const PEFFEntry& rhs) = default;
    PEFFEntry(PEFFEntry&& rhs) noexcept = default;
    PEFFEntry& operator=(const PEFFEntry& rhs) = default;
    PEFFEntry& operator=(PEFFEntry&& rhs) noexcept = default;

    bool operator==(const PEFFEntry& rhs) const
    {
      return prefix == rhs.prefix &&
             identifier == rhs.identifier &&
             sequence == rhs.sequence &&
             protein_names == rhs.protein_names &&
             gene_name == rhs.gene_name &&
             ncbi_tax_id == rhs.ncbi_tax_id &&
             taxonomy_name == rhs.taxonomy_name &&
             sequence_length == rhs.sequence_length &&
             sequence_version == rhs.sequence_version &&
             entry_version == rhs.entry_version &&
             protein_existence == rhs.protein_existence &&
             db_unique_id == rhs.db_unique_id &&
             entry_id == rhs.entry_id &&
             alt_accessions == rhs.alt_accessions &&
             modifications == rhs.modifications &&
             simple_variants == rhs.simple_variants &&
             complex_variants == rhs.complex_variants &&
             processed_regions == rhs.processed_regions &&
             disulfide_bonds == rhs.disulfide_bonds &&
             proteoforms == rhs.proteoforms &&
             custom_annotations == rhs.custom_annotations;
    }

    /// Convert to a FASTAFile::FASTAEntry (loses PEFF-specific annotations)
    FASTAFile::FASTAEntry toFASTAEntry() const;

    /// Create a PEFFEntry from a FASTAEntry (basic fields only)
    static PEFFEntry fromFASTAEntry(const FASTAFile::FASTAEntry& fasta);

    /**
      @brief Get the base AASequence for this entry (unmodified sequence).

      @return AASequence representing the protein sequence
    */
    AASequence getSequence() const;

    /**
      @brief Get an AASequence with all annotated modifications applied.

      Uses the modifications vector to apply modifications to the sequence.
      Modifications with unknown positions (position == 0) are skipped.
      Modifications that cannot be resolved are logged as warnings and skipped.

      @return AASequence with modifications applied
    */
    AASequence getModifiedSequence() const;

    /**
      @brief Get all variant sequences (each variant applied individually).

      @param descriptions Output vector for variant descriptions
      @param sequences Output vector for variant sequences
      @param include_complex If true, also include complex variants (default: false)
    */
    void getVariantSequences(
      std::vector<std::string>& descriptions,
      std::vector<AASequence>& sequences,
      bool include_complex = false) const;

    /**
      @brief Get processed sequence (e.g., mature protein without signal peptide).

      Applies the first processed region of the given type to extract
      the processed sequence segment.

      @param region_accession PEFF CV accession for the region type (e.g., "PEFF:0001021" for signal peptide)
      @return Processed AASequence, or empty if region not found
    */
    AASequence getProcessedSequence(const String& region_accession = "PEFF:0001021") const;

    /**
      @brief Generate all variant and/or modification peptides by digesting with a given protease.

      This method performs enzymatic digestion on the reference sequence and then
      generates all combinations of simple variants and/or modifications within each peptide.

      @param digestor The protease digestion object (must be configured with enzyme, missed cleavages, etc.)
      @param descriptions Output vector for peptide descriptions (empty for reference peptides)
      @param sequences Output vector for peptide sequences
      @param min_length Minimum peptide length to include (default: 6)
      @param max_length Maximum peptide length to include (default: 40, 0 = no limit)
      @param include_reference If true, include reference peptides (default: true)
      @param include_variants If true, generate variant combinations (default: true)
      @param include_modifications If true, generate modification combinations (default: false)
    */
    void digestWithVariants(
      const ProteaseDigestion& digestor,
      std::vector<std::string>& descriptions,
      std::vector<AASequence>& sequences,
      Size min_length = 6,
      Size max_length = 40,
      bool include_reference = true,
      bool include_variants = true,
      bool include_modifications = false) const;

    /**
      @brief Generate peptides with PEFF annotations and optional sample handling modifications.

      Combines enzymatic digestion, PEFF variants/modifications, and sample handling mods.

      @param digestor The protease digestion object (configured with enzyme, missed cleavages)
      @param descriptions Output vector for peptide descriptions
      @param sequences Output vector for peptide sequences
      @param fixed_mods Fixed modifications (e.g., {"Carbamidomethyl (C)"})
      @param variable_mods Variable modifications (e.g., {"Oxidation (M)"})
      @param max_variable_mods_per_peptide Maximum variable mods per peptide (default: 2)
      @param min_length Minimum peptide length (default: 6)
      @param max_length Maximum peptide length (default: 40, 0 = no limit)
      @param include_reference Include reference peptides (default: true)
      @param include_peff_variants Enumerate PEFF variants (default: true)
      @param include_peff_modifications Enumerate PEFF modifications (default: true)
    */
    void generatePeptides(
      const ProteaseDigestion& digestor,
      std::vector<std::string>& descriptions,
      std::vector<AASequence>& sequences,
      const std::vector<std::string>& fixed_mods = {},
      const std::vector<std::string>& variable_mods = {},
      Size max_variable_mods_per_peptide = 2,
      Size min_length = 6,
      Size max_length = 40,
      bool include_reference = true,
      bool include_peff_variants = true,
      bool include_peff_modifications = true) const;

  private:
    /**
      @brief Apply PEFF modifications at specific positions to a peptide.

      Helper method that generates all 2^n combinations of PEFF modifications
      for a given peptide, where n is the number of PEFF modifications within
      the peptide's range.

      @param peptide The base peptide sequence
      @param peff_mods PEFF modifications with positions relative to the peptide (0-based)
      @param base_description Base description to prepend to modification descriptions
      @return Vector of pairs: (description, modified AASequence)
    */
    static std::vector<std::pair<String, AASequence>> enumeratePEFFModifications_(
      const AASequence& peptide,
      const std::vector<std::pair<Size, const PEFFModification*>>& peff_mods,
      const String& base_description);
  };

  /**
    @brief Represents a custom key definition from the PEFF header.
  */
  struct OPENMS_DLLAPI PEFFCustomKeyDef
  {
    String key_name;
    String description;
    String concept_curie;
    String regexp;
    std::vector<String> field_names;
    std::vector<String> field_types;

    bool operator==(const PEFFCustomKeyDef& rhs) const
    {
      return key_name == rhs.key_name && description == rhs.description &&
             concept_curie == rhs.concept_curie && regexp == rhs.regexp &&
             field_names == rhs.field_names && field_types == rhs.field_types;
    }
  };

  /**
    @brief Metadata from a PEFF database header section.

    The header section contains lines starting with # that describe the database.
  */
  struct OPENMS_DLLAPI PEFFDatabaseMetadata
  {
    String version{"1.0"};
    String db_name;
    String prefix;
    String db_description;
    bool is_decoy{false};
    std::vector<String> db_sources;
    String db_version;
    Size number_of_entries{0};

    enum class SequenceType { AA, NA };
    SequenceType sequence_type{SequenceType::AA};

    std::vector<String> general_comments;     ///< Multiple GeneralComment lines allowed
    String conversion;                         ///< Conversion notes
    std::map<String, String> specific_keys;    ///< SpecificKey definitions (key -> description)
    std::map<String, String> specific_values;  ///< SpecificValue definitions (key -> type)
    std::vector<String> optional_tag_defs;
    std::vector<PEFFCustomKeyDef> custom_key_defs;  ///< CustomKeyDef definitions
    bool has_annotation_identifiers{false};    ///< Whether entries use annotation identifiers
    bool is_proteoform_db{false};              ///< Whether this is a proteoform database (ProteoformDb)

    std::map<String, String> unrecognized_keys;  ///< Unrecognized header keys (preserved for round-trip)

    PEFFDatabaseMetadata() = default;

    bool operator==(const PEFFDatabaseMetadata& rhs) const
    {
      return version == rhs.version &&
             db_name == rhs.db_name &&
             prefix == rhs.prefix &&
             db_description == rhs.db_description &&
             is_decoy == rhs.is_decoy &&
             db_sources == rhs.db_sources &&
             db_version == rhs.db_version &&
             number_of_entries == rhs.number_of_entries &&
             sequence_type == rhs.sequence_type &&
             general_comments == rhs.general_comments &&
             conversion == rhs.conversion &&
             specific_keys == rhs.specific_keys &&
             specific_values == rhs.specific_values &&
             optional_tag_defs == rhs.optional_tag_defs &&
             custom_key_defs == rhs.custom_key_defs &&
             has_annotation_identifiers == rhs.has_annotation_identifiers &&
             is_proteoform_db == rhs.is_proteoform_db &&
             unrecognized_keys == rhs.unrecognized_keys;
    }
  };

  /**
    @brief This class serves for reading and writing PEFF (PSI Extended FASTA Format) files.

    PEFF extends FASTA with rich annotations for modifications, variants, processed regions,
    and proteoforms. See https://github.com/HUPO-PSI/PEFF for the specification.

    You can use aggregate methods load() and store() to read/write a set of protein
    sequences at the cost of memory.

    Or use single read/write of protein sequences using readStart(), readNext()
    and writeStart(), writeNext(), writeEnd() for more memory efficiency.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI PEFFFile : public ProgressLogger
  {
  public:
    /// Default constructor
    PEFFFile() = default;

    /// Destructor
    ~PEFFFile() override = default;

    /**
      @brief Loads a PEFF file and stores entries and headers.

      @param filename The PEFF file to load
      @param entries Output vector for PEFF entries
      @param headers Output vector for database metadata (one per database in file)

      @exception Exception::FileNotFound is thrown if the file does not exist.
      @exception Exception::ParseError is thrown if the file format is invalid.
    */
    void load(const String& filename,
              std::vector<PEFFEntry>& entries,
              std::vector<PEFFDatabaseMetadata>& headers) const;

    /**
      @brief Stores entries to a PEFF file with the given header.

      @param filename The output file path
      @param entries The entries to store
      @param header The database metadata header

      @exception Exception::UnableToCreateFile is thrown if the file cannot be created.
    */
    void store(const String& filename,
               const std::vector<PEFFEntry>& entries,
               const PEFFDatabaseMetadata& header) const;

    /**
      @brief Stores entries to a PEFF file with multiple database headers.

      Writes a single file description block followed by one sequence database
      description block per provided header, then all entries.

      @param filename The output file path
      @param entries The entries to store
      @param headers The database metadata headers (one per database in file)

      @exception Exception::UnableToCreateFile is thrown if the file cannot be created.
      @exception Exception::InvalidParameter is thrown if headers is empty.
    */
    void store(const String& filename,
               const std::vector<PEFFEntry>& entries,
               const std::vector<PEFFDatabaseMetadata>& headers) const;

    /**
      @brief Prepares a PEFF file for streamed reading using readNext().

      @param filename The PEFF file to read

      @exception Exception::FileNotFound is thrown if the file does not exist.
      @exception Exception::FileNotReadable is thrown if the file cannot be read.
    */
    void readStart(const String& filename);

    /**
      @brief Reads the next PEFF entry from the file.

      @param entry Output for the next entry
      @return true if an entry was read, false if EOF was reached

      @exception Exception::ParseError is thrown if parsing fails.
    */
    bool readNext(PEFFEntry& entry);

    /**
      @brief Returns the headers parsed during readStart().

      Headers are available after calling readStart().
    */
    const std::vector<PEFFDatabaseMetadata>& getHeaders() const;

    /// Returns true if the end of the file has been reached
    bool atEnd() const;

    /**
      @brief Prepares a PEFF file for streamed writing using writeNext().

      @param filename The output file path
      @param header The database metadata header

      @exception Exception::UnableToCreateFile is thrown if the file cannot be created.
    */
    void writeStart(const String& filename, const PEFFDatabaseMetadata& header);

    /**
      @brief Prepares a PEFF file for streamed writing using writeNext(), with multiple headers.

      Writes a single file description block followed by one database description block
      per header.

      @param filename The output file path
      @param headers The database metadata headers (one per database in file)

      @exception Exception::UnableToCreateFile is thrown if the file cannot be created.
      @exception Exception::InvalidParameter is thrown if headers is empty.
    */
    void writeStart(const String& filename, const std::vector<PEFFDatabaseMetadata>& headers);

    /**
      @brief Writes the next PEFF entry to the file.

      @param entry The entry to write
    */
    void writeNext(const PEFFEntry& entry);

    /// Closes the output file (called automatically in destructor)
    void writeEnd();

    /**
      @brief Checks if a file appears to be a PEFF file (by checking for # PEFF header).

      @param filename The file to check
      @return true if the file starts with PEFF headers
    */
    static bool isPEFFFile(const String& filename);

    /**
      @brief Converts a PEFF entry to ProForma notation.

      @param entry The PEFF entry to convert
      @return ProForma string representation
    */
    static String toProForma(const PEFFEntry& entry);

  protected:
    /// Parse a header line (# Key=Value or # //)
    void parseHeaderLine_(const String& line, PEFFDatabaseMetadata& header, bool& new_db);

    /// Parse annotations from the description line
    void parseAnnotations_(const String& description, PEFFEntry& entry);

    /// Parse a single modification tuple
    PEFFModification parseModification_(const String& tuple);

    /// Parse a simple variant tuple
    PEFFVariantSimple parseVariantSimple_(const String& tuple);

    /// Parse a complex variant tuple
    PEFFVariantComplex parseVariantComplex_(const String& tuple);

    /// Parse a processed region tuple
    PEFFProcessedRegion parseProcessedRegion_(const String& tuple);

    /// Parse a disulfide bond tuple
    PEFFDisulfideBond parseDisulfideBond_(const String& tuple);

    /// Parse a parenthesized list of values
    std::vector<String> parseParenList_(const String& value);

    /// Format the header section for output
    String formatHeader_(const PEFFDatabaseMetadata& header) const;

    /// Format the header section for output (multiple database blocks)
    String formatHeader_(const std::vector<PEFFDatabaseMetadata>& headers) const;

    /// Format a single entry for output
    String formatEntry_(const PEFFEntry& entry) const;


    /// Read entry data (identifier, description, sequence)
    bool readEntry_(std::string& id, std::string& description, std::string& seq);

    std::fstream infile_;                           ///< Input file stream
    std::ofstream outfile_;                         ///< Output file stream
    std::vector<PEFFDatabaseMetadata> headers_;     ///< Parsed headers
    Size entries_read_{0};                          ///< Number of entries read
    std::streampos fileSize_{0};                    ///< File size for progress
    std::string seq_;                               ///< Current sequence buffer
    std::string id_;                                ///< Current identifier buffer
    std::string description_;                       ///< Current description buffer
  };

} // namespace OpenMS
