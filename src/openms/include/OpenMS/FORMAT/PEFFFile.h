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

#include <fstream>
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
    String evidence;          ///< Optional evidence tag
    String annotation_id;     ///< Optional annotation identifier (when HasAnnotationIdentifiers=true)

    enum class Type { PSI_MOD, UNIMOD, GENERIC };
    Type type{Type::GENERIC};

    PEFFModification() = default;
    PEFFModification(Size pos, const String& acc, const String& n, const String& ev = "", const String& aid = "")
      : position(pos), accession(acc), name(n), evidence(ev), annotation_id(aid)
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
             name == rhs.name && evidence == rhs.evidence && type == rhs.type &&
             annotation_id == rhs.annotation_id;
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
    String sources;           ///< Source references (dbSNP, COSMIC, etc.)
    String annotation_id;     ///< Optional annotation identifier (when HasAnnotationIdentifiers=true)

    PEFFVariantSimple() = default;
    PEFFVariantSimple(Size pos, char aa, const String& src = "", const String& aid = "")
      : position(pos), variant_aa(aa), sources(src), annotation_id(aid) {}

    bool operator==(const PEFFVariantSimple& rhs) const
    {
      return position == rhs.position && variant_aa == rhs.variant_aa &&
             sources == rhs.sources && annotation_id == rhs.annotation_id;
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
    String sources;           ///< Source references
    String annotation_id;     ///< Optional annotation identifier (when HasAnnotationIdentifiers=true)

    PEFFVariantComplex() = default;
    PEFFVariantComplex(Size start, Size end, const String& repl, const String& src = "", const String& aid = "")
      : start_position(start), end_position(end), replacement(repl), sources(src), annotation_id(aid) {}

    bool operator==(const PEFFVariantComplex& rhs) const
    {
      return start_position == rhs.start_position && end_position == rhs.end_position &&
             replacement == rhs.replacement && sources == rhs.sources &&
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
    String type;              ///< PEFF CV term (e.g., "PEFF:0001021")
    String name;              ///< Optional name (e.g., "signal peptide")
    String description;       ///< Optional description
    String annotation_id;     ///< Optional annotation identifier (when HasAnnotationIdentifiers=true)

    PEFFProcessedRegion() = default;
    PEFFProcessedRegion(Size start, Size end, const String& t, const String& n = "", const String& desc = "", const String& aid = "")
      : start_position(start), end_position(end), type(t), name(n), description(desc), annotation_id(aid) {}

    bool operator==(const PEFFProcessedRegion& rhs) const
    {
      return start_position == rhs.start_position && end_position == rhs.end_position &&
             type == rhs.type && name == rhs.name && description == rhs.description &&
             annotation_id == rhs.annotation_id;
    }
  };

  /**
    @brief Represents a single entry in a PEFF file with all annotations.
  */
  struct OPENMS_DLLAPI PEFFEntry
  {
    // Basic fields
    String identifier;
    String sequence;

    // Metadata
    std::vector<String> protein_names;   ///< \\PName - may have multiple names
    String gene_name;                     ///< \\GName
    Int ncbi_tax_id{0};                  ///< \\NcbiTaxId
    String taxonomy_name;                 ///< \\TaxName
    Size sequence_length{0};             ///< \\Length
    String sequence_version;              ///< \\SV
    String entry_version;                 ///< \\EV
    Int protein_existence{0};            ///< \\PE (1-5)

    // Annotations
    std::vector<PEFFModification> modifications;
    std::vector<PEFFVariantSimple> simple_variants;
    std::vector<PEFFVariantComplex> complex_variants;
    std::vector<PEFFProcessedRegion> processed_regions;
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
      return identifier == rhs.identifier &&
             sequence == rhs.sequence &&
             protein_names == rhs.protein_names &&
             gene_name == rhs.gene_name &&
             ncbi_tax_id == rhs.ncbi_tax_id &&
             taxonomy_name == rhs.taxonomy_name &&
             sequence_length == rhs.sequence_length &&
             sequence_version == rhs.sequence_version &&
             entry_version == rhs.entry_version &&
             protein_existence == rhs.protein_existence &&
             modifications == rhs.modifications &&
             simple_variants == rhs.simple_variants &&
             complex_variants == rhs.complex_variants &&
             processed_regions == rhs.processed_regions &&
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

      Each variant sequence has one variant applied.
      Does not combine variants.

      @param include_complex If true, also include complex variants (insertions, deletions, multi-residue substitutions)
      @return Vector of pairs: (variant description, AASequence)
    */
    std::vector<std::pair<String, AASequence>> getVariantSequences(bool include_complex = false) const;

    /**
      @brief Get processed sequence (e.g., mature protein without signal peptide).

      Applies the first processed region of the given type to extract
      the processed sequence segment.

      @param region_type PEFF CV term for the region type (e.g., "PEFF:0001021" for signal peptide)
      @return Processed AASequence, or empty if region not found
    */
    AASequence getProcessedSequence(const String& region_type = "PEFF:0001021") const;
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
    String db_date;                      ///< YYYYMMDD format
    Size number_of_entries{0};

    enum class SequenceType { AA, NA };
    SequenceType sequence_type{SequenceType::AA};

    String general_comment;
    std::map<String, String> specific_keys;
    std::vector<String> optional_tag_defs;
    bool has_annotation_identifiers{false};  ///< Whether entries use annotation identifiers

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
             db_date == rhs.db_date &&
             number_of_entries == rhs.number_of_entries &&
             sequence_type == rhs.sequence_type &&
             general_comment == rhs.general_comment &&
             specific_keys == rhs.specific_keys &&
             optional_tag_defs == rhs.optional_tag_defs &&
             has_annotation_identifiers == rhs.has_annotation_identifiers;
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
      @brief Prepares a PEFF file for streamed reading using readNext().

      @param filename The PEFF file to read

      @exception Exception::FileNotFound is thrown if the file does not exist.
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

    /// Parse a parenthesized list of values
    std::vector<String> parseParenList_(const String& value);

    /// Format the header section for output
    String formatHeader_(const PEFFDatabaseMetadata& header) const;

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
