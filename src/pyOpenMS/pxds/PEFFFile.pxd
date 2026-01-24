from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp.pair cimport pair as libcpp_pair
from libcpp cimport bool

from String cimport *
from FASTAFile cimport *
from AASequence cimport *

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


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFModification:
        # wrap-doc:
        #  Represents a PEFF modification annotation
        #
        #  Attributes:
        #    position: 1-based position (0 = unknown)
        #    accession: Modification accession (MOD:xxxxx, UNIMOD:xx, or custom)
        #    name: Human-readable name
        #    evidence: Optional evidence tag

        PEFFModification() except + nogil
        PEFFModification(PEFFModification) except + nogil

        Size position
        String accession
        String name
        String evidence


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFVariantSimple:
        # wrap-doc:
        #  Represents a simple PEFF variant (single amino acid substitution)
        #
        #  Attributes:
        #    position: 1-based position
        #    variant_aa: Variant amino acid character
        #    sources: Source references (dbSNP, COSMIC, etc.)

        PEFFVariantSimple() except + nogil
        PEFFVariantSimple(PEFFVariantSimple) except + nogil

        Size position
        char variant_aa
        String sources


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFVariantComplex:
        # wrap-doc:
        #  Represents a complex PEFF variant (multi-residue change)
        #
        #  Attributes:
        #    start_position: 1-based start position
        #    end_position: 1-based end position
        #    replacement: Replacement sequence (empty = deletion)
        #    sources: Source references

        PEFFVariantComplex() except + nogil
        PEFFVariantComplex(PEFFVariantComplex) except + nogil

        Size start_position
        Size end_position
        String replacement
        String sources


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFProcessedRegion:
        # wrap-doc:
        #  Represents a PEFF processed region (signal peptide, etc.)
        #
        #  Attributes:
        #    start_position: 1-based start position
        #    end_position: 1-based end position
        #    type: PEFF CV term (e.g., PEFF:0001021)
        #    name: Optional name
        #    description: Optional description

        PEFFProcessedRegion() except + nogil
        PEFFProcessedRegion(PEFFProcessedRegion) except + nogil

        Size start_position
        Size end_position
        String type
        String name
        String description


cdef extern from "<OpenMS/FORMAT/PEFFFile.h>" namespace "OpenMS":

    cdef cppclass PEFFEntry:
        # wrap-doc:
        #  Represents a single entry in a PEFF file with all annotations
        #
        #  Attributes:
        #    identifier: Protein identifier
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

        libcpp_vector[PEFFModification] modifications
        libcpp_vector[PEFFVariantSimple] simple_variants
        libcpp_vector[PEFFVariantComplex] complex_variants
        libcpp_vector[PEFFProcessedRegion] processed_regions
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

        libcpp_vector[libcpp_pair[String, AASequence]] getVariantSequences() except + nogil
            # wrap-doc:
            #  Get all variant sequences (each simple variant applied individually).
            #  Each variant sequence has one amino acid substitution applied.
            #  Returns vector of pairs: (variant description, AASequence)

        AASequence getProcessedSequence(String region_type) except + nogil
            # wrap-doc:
            #  Get processed sequence (e.g., mature protein without signal peptide).
            #  For signal peptides (PEFF:0001021), returns the mature protein.
            #  :param region_type: PEFF CV term (default "PEFF:0001021" for signal peptide)


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
        #    db_date: Database date (YYYYMMDD)
        #    number_of_entries: Number of entries
        #    sequence_type: Sequence type (AA or NA)
        #    general_comment: General comment

        PEFFDatabaseMetadata() except + nogil
        PEFFDatabaseMetadata(PEFFDatabaseMetadata) except + nogil

        String version
        String db_name
        String prefix
        String db_description
        bool is_decoy
        libcpp_vector[String] db_sources
        String db_version
        String db_date
        Size number_of_entries
        SequenceType sequence_type
        String general_comment
        libcpp_map[String, String] specific_keys
        libcpp_vector[String] optional_tag_defs
