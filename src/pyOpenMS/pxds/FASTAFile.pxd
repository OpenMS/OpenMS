from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool

from String cimport *
from DefaultParamHandler cimport *

from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/FORMAT/FASTAFile.h>" namespace "OpenMS":

    cdef cppclass FASTAFile:
        # wrap-doc:
        #  File adapter for FASTA files
        #  
        #  Provides methods to load and store protein/peptide sequences in FASTA format.
        #  Supports both batch loading (load/store) and streaming (readStart/readNext/writeStart/writeNext).
        #  
        #  Usage:
        #  
        #  .. code-block:: python
        #  
        #    # Batch loading
        #    entries = []
        #    FASTAFile().load("proteins.fasta", entries)
        #    for entry in entries:
        #      print(entry.identifier, entry.sequence)
        #  
        #    # Streaming (memory-efficient for large files)
        #    ff = FASTAFile()
        #    ff.readStart("proteins.fasta")
        #    entry = FASTAEntry()
        #    while ff.readNext(entry):
        #      print(entry.identifier)

        FASTAFile() except + nogil
        # copy constructor of 'FASTAFile' is implicitly deleted because field 'infile_' has a deleted copy constructor
        FASTAFile(FASTAFile &) except + nogil  # wrap-ignore

        void load(const String& filename, libcpp_vector[FASTAEntry] & data) except + nogil  # wrap-doc:Loads a FASTA file given by 'filename' and stores the information in 'data'
        void store(const String& filename, libcpp_vector[FASTAEntry] & data) except + nogil  # wrap-doc:Stores the data given by 'data' at the file 'filename'

        void readStart(const String & filename) except + nogil
            # wrap-doc:
            #  Prepares a FASTA file given by 'filename' for streamed reading using readNext()
            #  
            #  :raises:
            #      Exception:FileNotFound is thrown if the file does not exists
            #  :raises:
            #      Exception:ParseError is thrown if the file does not suit to the standard
        bool readNext(FASTAEntry & protein) except + nogil
            # wrap-doc:
            #  Reads the next FASTA entry from file
            #  
            #  If you want to read all entries in one go, use load()
            #  
            #  :return: true if entry was read; false if eof was reached
            #  :raises:
            #      Exception:FileNotFound is thrown if the file does not exists
            #  :raises:
            #      Exception:ParseError is thrown if the file does not suit to the standard
        # NAMESPACE # std::streampos position() except + nogil
        bool atEnd() except + nogil  # wrap-doc:Boolean function to check if streams is at end of file
        # NAMESPACE # bool setPosition(const std::streampos & pos) except + nogil
        void writeStart(const String & filename) except + nogil
            # wrap-doc:
            #  Prepares a FASTA file given by 'filename' for streamed writing using writeNext()
            #  
            #  :raises:
            #      Exception:UnableToCreateFile is thrown if the process is not able to write to the file (disk full?)
        void writeNext(const FASTAEntry & protein) except + nogil
            # wrap-doc:
            #  Stores the data given by `protein`. Call writeStart() once before calling writeNext()
            #  
            #  Call writeEnd() when done to close the file!
            #  
            #  :raises:
            #      Exception:UnableToCreateFile is thrown if the process is not able to write to the file (disk full?)
        void writeEnd() except + nogil  # wrap-doc:Closes the file (flush). Called implicitly when FASTAFile object does out of scope

cdef extern from "<OpenMS/FORMAT/FASTAFile.h>" namespace "OpenMS::FASTAFile":

    cdef cppclass FASTAEntry:
        # wrap-doc:
        #  Represents a single FASTA entry with identifier, description, and sequence
        #  
        #  Attributes:
        #    identifier: The protein/sequence identifier (from FASTA header before first space)
        #    description: The description line (from FASTA header after first space)
        #    sequence: The amino acid or nucleotide sequence

        FASTAEntry() except + nogil
        FASTAEntry(FASTAEntry) except + nogil

        String identifier  # wrap-doc:Protein/sequence identifier
        String description  # wrap-doc:Description from FASTA header
        String sequence  # wrap-doc:The amino acid or nucleotide sequence

        bool headerMatches(const FASTAEntry & rhs) except + nogil  # wrap-doc:Returns True if identifier and description match the given entry
        bool sequenceMatches(const FASTAEntry & rhs) except + nogil  # wrap-doc:Returns True if the sequence matches the given entry
