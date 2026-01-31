from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from libcpp.string cimport string as libcpp_utf8_string
from libcpp.string cimport string as libcpp_utf8_output_string
from libcpp cimport bool
from Types cimport *

from String cimport *
from FASTAFile cimport *
from AASequence cimport *
from ProteaseDigestion cimport *

cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFFile:
        # wrap-doc:
        #  File adapter for PEFF (PSI Extended FASTA Format) files
        #
        #  PEFF extends FASTA with rich annotations for modifications, variants,
        #  processed regions, and proteoforms.
        #
        #  Usage:
        #
        #  .. code-block:: python
        #
        #    # Batch loading
        #    peff = PEFFFile()
        #    entries = []
        #    headers = []
        #    peff.load("proteins.peff", entries, headers)
        #    for entry in entries:
        #      print(entry.identifier, len(entry.modifications))
        #
        #    # Streaming (memory-efficient for large files)
        #    peff = PEFFFile()
        #    peff.readStart("proteins.peff")
        #    entry = PEFFEntry()
        #    while peff.readNext(entry):
        #      print(entry.identifier)

        PEFFFile() except + nogil
        PEFFFile(PEFFFile &) except + nogil  # wrap-ignore

        void load(const String& filename,
                  libcpp_vector[PEFFEntry] & entries,
                  libcpp_vector[PEFFDatabaseMetadata] & headers) except + nogil
            # wrap-doc:
            #  Loads a PEFF file and stores entries and headers
            #
            #  :param filename: The PEFF file to load
            #  :param entries: Output list for PEFF entries
            #  :param headers: Output list for database metadata

        void store(const String& filename,
                   const libcpp_vector[PEFFEntry] & entries,
                   const PEFFDatabaseMetadata & header) except + nogil
            # wrap-doc:
            #  Stores entries to a PEFF file with the given header
            #
            #  :param filename: The output file path
            #  :param entries: The entries to store
            #  :param header: The database metadata header

        void readStart(const String & filename) except + nogil
            # wrap-doc:
            #  Prepares a PEFF file for streamed reading using readNext()
            #
            #  :raises:
            #      Exception:FileNotFound is thrown if the file does not exist
            #      Exception:FileNotReadable is thrown if the file cannot be read

        bool readNext(PEFFEntry & entry) except + nogil
            # wrap-doc:
            #  Reads the next PEFF entry from file
            #
            #  :return: True if entry was read; False if EOF was reached

        libcpp_vector[PEFFDatabaseMetadata] getHeaders() except + nogil
            # wrap-doc:
            #  Returns the headers parsed during readStart()

        bool atEnd() except + nogil
            # wrap-doc:
            #  Returns True if the end of the file has been reached

        void writeStart(const String & filename, const PEFFDatabaseMetadata & header) except + nogil
            # wrap-doc:
            #  Prepares a PEFF file for streamed writing using writeNext()
            #
            #  :raises:
            #      Exception:UnableToCreateFile is thrown if the file cannot be created

        void writeNext(const PEFFEntry & entry) except + nogil
            # wrap-doc:
            #  Writes the next PEFF entry to the file

        void writeEnd() except + nogil
            # wrap-doc:
            #  Closes the output file

    bool isPEFFFile "OpenMS::PEFFFile::isPEFFFile" (const String& filename) except + nogil
        # wrap-doc:
        #  Checks if a file appears to be a PEFF file
        #
        #  :param filename: The file to check
        #  :return: True if the file starts with PEFF headers

    String toProForma "OpenMS::PEFFFile::toProForma" (const PEFFEntry& entry) except + nogil
        # wrap-doc:
        #  Converts a PEFF entry to ProForma notation
        #
        #  :param entry: The PEFF entry to convert
        #  :return: ProForma string representation


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS::PEFFModification":

    cdef enum class ModificationType "OpenMS::PEFFModification::Type":
        # wrap-attach:
        #    PEFFModification
        PSI_MOD
        UNIMOD
        GENERIC


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFModification:
        # wrap-doc:
        #  Represents a PEFF modification annotation
        #
        #  Attributes:
        #    position: 1-based position (0 = unknown)
        #    accession: Modification accession (MOD:xxxxx, UNIMOD:xx, or custom)
        #    name: Human-readable name
        #    optional_tag: Optional tag (last component of annotation tuple)
        #    annotation_id: Optional annotation identifier (UInt, max value = not set)
        #    type: Modification type (PSI_MOD, UNIMOD, or GENERIC)

        PEFFModification() except + nogil
        PEFFModification(PEFFModification) except + nogil

        Size position
        String accession
        String name
        String optional_tag
        UInt annotation_id
        ModificationType type


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFVariantSimple:
        # wrap-doc:
        #  Represents a simple PEFF variant (single amino acid substitution)
        #
        #  Attributes:
        #    position: 1-based position
        #    variant_aa: Variant amino acid character
        #    optional_tag: Optional tag (last component of annotation tuple)
        #    annotation_id: Optional annotation identifier (UInt, max value = not set)

        PEFFVariantSimple() except + nogil
        PEFFVariantSimple(PEFFVariantSimple) except + nogil

        Size position
        char variant_aa
        String optional_tag
        UInt annotation_id


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFVariantComplex:
        # wrap-doc:
        #  Represents a complex PEFF variant (multi-residue change)
        #
        #  Attributes:
        #    start_position: 1-based start position
        #    end_position: 1-based end position
        #    replacement: Replacement sequence (empty = deletion)
        #    optional_tag: Optional tag (last component of annotation tuple)
        #    annotation_id: Optional annotation identifier (UInt, max value = not set)

        PEFFVariantComplex() except + nogil
        PEFFVariantComplex(PEFFVariantComplex) except + nogil

        Size start_position
        Size end_position
        String replacement
        String optional_tag
        UInt annotation_id


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFProcessedRegion:
        # wrap-doc:
        #  Represents a PEFF processed region (signal peptide, etc.)
        #
        #  Attributes:
        #    start_position: 1-based start position
        #    end_position: 1-based end position
        #    accession: PEFF CV accession (e.g., PEFF:0001021)
        #    name: Optional name
        #    optional_tag: Optional tag (last component of annotation tuple)
        #    annotation_id: Optional annotation identifier (UInt, max value = not set)

        PEFFProcessedRegion() except + nogil
        PEFFProcessedRegion(PEFFProcessedRegion) except + nogil

        Size start_position
        Size end_position
        String accession
        String name
        String optional_tag
        UInt annotation_id


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFDisulfideBond:
        # wrap-doc:
        #  Represents a disulfide bond annotation in PEFF
        #
        #  Attributes:
        #    id1: First cysteine reference (AnnotationIdentifier)
        #    id2: Second cysteine reference (AnnotationIdentifier)
        #    optional_tag: Optional tag (e.g., "between chains")
        #    annotation_id: Optional annotation identifier (UInt, max value = not set)

        PEFFDisulfideBond() except + nogil
        PEFFDisulfideBond(PEFFDisulfideBond) except + nogil

        String id1
        String id2
        String optional_tag
        UInt annotation_id


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFCustomKeyDef:
        # wrap-doc:
        #  Represents a custom key definition from the PEFF header
        #
        #  Attributes:
        #    key_name: Custom key name
        #    description: Description of the key

        PEFFCustomKeyDef() except + nogil
        PEFFCustomKeyDef(PEFFCustomKeyDef) except + nogil

        String key_name
        String description
        String concept_curie
        String regexp
        libcpp_vector[String] field_names
        libcpp_vector[String] field_types


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFEntry:
        # wrap-doc:
        #  Represents a single entry in a PEFF file with all annotations.
        #
        #  Each entry corresponds to one description line and sequence.
        #  The description line format is: >Prefix:DbUniqueId \\key=value ...
        #
        #  Attributes:
        #    prefix: Database prefix from description line
        #    identifier: Protein identifier (Prefix:DbUniqueId)
        #    sequence: Amino acid sequence
        #    protein_names: List of protein names
        #    gene_name: Gene name
        #    ncbi_tax_id: NCBI taxonomy ID
        #    taxonomy_name: Taxonomy name
        #    sequence_length: Sequence length
        #    modifications: List of modifications
        #    simple_variants: List of simple variants
        #    complex_variants: List of complex variants
        #    processed_regions: List of processed regions
        #    proteoforms: List of proteoforms in ProForma notation

        PEFFEntry() except + nogil
        PEFFEntry(PEFFEntry) except + nogil

        String prefix
        String identifier
        String sequence
        libcpp_vector[String] protein_names
        String gene_name
        int ncbi_tax_id
        String taxonomy_name
        Size sequence_length
        String sequence_version
        String entry_version
        int protein_existence
        String db_unique_id
        String entry_id
        libcpp_vector[String] alt_accessions

        libcpp_vector[PEFFModification] modifications
        libcpp_vector[PEFFVariantSimple] simple_variants
        libcpp_vector[PEFFVariantComplex] complex_variants
        libcpp_vector[PEFFProcessedRegion] processed_regions
        libcpp_vector[PEFFDisulfideBond] disulfide_bonds
        libcpp_vector[String] proteoforms
        libcpp_map[String, String] custom_annotations

        FASTAEntry toFASTAEntry() except + nogil
            # wrap-doc:
            #  Convert to a FASTAEntry (loses PEFF-specific annotations)

        AASequence getSequence() except + nogil
            # wrap-doc:
            #  Get the base AASequence for this entry (unmodified sequence)

        AASequence getModifiedSequence() except + nogil
            # wrap-doc:
            #  Get an AASequence with all annotated modifications applied.
            #  Modifications with unknown positions (position == 0) are skipped.

        void getVariantSequences(
            libcpp_vector[libcpp_utf8_output_string]& descriptions,
            libcpp_vector[AASequence]& sequences,
            bool include_complex) except + nogil
            # wrap-doc:
            #  Get all variant sequences (each variant applied individually).
            #  :param descriptions: Output list for variant descriptions
            #  :param sequences: Output list for variant sequences
            #  :param include_complex: Include complex variants (default: false)

        AASequence getProcessedSequence() except + nogil
            # wrap-doc:
            #  Get processed sequence (e.g., mature protein without signal peptide).
            #  Uses default region accession "PEFF:0001021" (signal peptide).

        AASequence getProcessedSequence(String region_accession) except + nogil
            # wrap-doc:
            #  Get processed sequence for the given region accession.
            #  :param region_accession: PEFF CV accession (e.g., "PEFF:0001021" for signal peptide)

        void digestWithVariants(
            const ProteaseDigestion& digestor,
            libcpp_vector[libcpp_utf8_output_string]& descriptions,
            libcpp_vector[AASequence]& sequences,
            Size min_length,
            Size max_length,
            bool include_reference,
            bool include_variants,
            bool include_modifications) except + nogil
            # wrap-doc:
            #  Digest protein and enumerate PEFF variant/modification combinations.
            #  :param digestor: ProteaseDigestion object
            #  :param descriptions: Output list for peptide descriptions
            #  :param sequences: Output list for peptide sequences

        void generatePeptides(
            const ProteaseDigestion& digestor,
            libcpp_vector[libcpp_utf8_output_string]& descriptions,
            libcpp_vector[AASequence]& sequences,
            const libcpp_vector[libcpp_utf8_string]& fixed_mods,
            const libcpp_vector[libcpp_utf8_string]& variable_mods,
            Size max_variable_mods_per_peptide,
            Size min_length,
            Size max_length,
            bool include_reference,
            bool include_peff_variants,
            bool include_peff_modifications) except + nogil
            # wrap-doc:
            #  Generate peptides with PEFF variants and sample modifications.
            #  :param digestor: ProteaseDigestion object
            #  :param descriptions: Output list for peptide descriptions
            #  :param sequences: Output list for peptide sequences
            #  :param fixed_mods: Fixed mods like ["Carbamidomethyl (C)"]
            #  :param variable_mods: Variable mods like ["Oxidation (M)"]


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS::PEFFEntry":

    PEFFEntry fromFASTAEntry "OpenMS::PEFFEntry::fromFASTAEntry" (const FASTAEntry&) except + nogil
        # wrap-doc:
        #  Create a PEFFEntry from a FASTAEntry


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS::PEFFDatabaseMetadata":

    cdef enum class SequenceType "OpenMS::PEFFDatabaseMetadata::SequenceType":
        # wrap-attach:
        #    PEFFDatabaseMetadata
        AA
        NA


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFDatabaseMetadata:
        # wrap-doc:
        #  Metadata from a PEFF database header section
        #
        #  Attributes:
        #    version: PEFF version
        #    db_name: Database name
        #    prefix: Database prefix
        #    db_description: Database description
        #    is_decoy: Whether this is a decoy database
        #    db_sources: List of database sources
        #    db_version: Database version
        #    number_of_entries: Number of entries
        #    sequence_type: Sequence type (AA or NA)
        #    general_comments: List of general comments
        #    conversion: Conversion notes
        #    specific_keys: Custom key definitions
        #    specific_values: Custom value type definitions
        #    has_annotation_identifiers: Whether entries use annotation identifiers
        #    is_proteoform_db: Whether this is a proteoform database
        #    unrecognized_keys: Unrecognized header keys preserved for round-trip

        PEFFDatabaseMetadata() except + nogil
        PEFFDatabaseMetadata(PEFFDatabaseMetadata) except + nogil

        String version
        String db_name
        String prefix
        String db_description
        bool is_decoy
        libcpp_vector[String] db_sources
        String db_version
        Size number_of_entries
        SequenceType sequence_type
        libcpp_vector[String] general_comments
        String conversion
        libcpp_map[String, String] specific_keys
        libcpp_map[String, String] specific_values
        libcpp_vector[String] optional_tag_defs
        libcpp_vector[PEFFCustomKeyDef] custom_key_defs
        bool has_annotation_identifiers
        bool is_proteoform_db
        libcpp_map[String, String] unrecognized_keys
