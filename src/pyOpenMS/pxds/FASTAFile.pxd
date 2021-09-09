from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool

from String cimport *
from DefaultParamHandler cimport *

from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/FORMAT/FASTAFile.h>" namespace "OpenMS":

    cdef cppclass FASTAFile:

        FASTAFile() nogil except + # wrap-doc:This class serves for reading in and writing FASTA files
        # copy constructor of 'FASTAFile' is implicitly deleted because field 'infile_' has a deleted copy constructor
        FASTAFile(FASTAFile &) nogil except + # wrap-ignore

        void load(const String& filename, libcpp_vector[FASTAEntry] & data) nogil except + # wrap-doc:Loads a FASTA file given by 'filename' and stores the information in 'data'
        void store(const String& filename, libcpp_vector[FASTAEntry] & data) nogil except + # wrap-doc:Stores the data given by 'data' at the file 'filename'

        void readStart(const String & filename) nogil except + 
            # wrap-doc:
            #   Prepares a FASTA file given by 'filename' for streamed reading using readNext()
            #   -----
            #   :raises:
            #       Exception:FileNotFound is thrown if the file does not exists
            #   :raises:
            #       Exception:ParseError is thrown if the file does not suit to the standard
        bool readNext(FASTAEntry & protein) nogil except +
            # wrap-doc:
            #   Reads the next FASTA entry from file
            #   -----
            #   If you want to read all entries in one go, use load()
            #   -----
            #   :returns: true if entry was read; false if eof was reached
            #   :raises:
            #       Exception:FileNotFound is thrown if the file does not exists
            #   :raises:
            #       Exception:ParseError is thrown if the file does not suit to the standard
        # NAMESPACE # std::streampos position() nogil except +
        bool atEnd() nogil except + # wrap-doc:Boolean function to check if streams is at end of file
        # NAMESPACE # bool setPosition(const std::streampos & pos) nogil except +
        void writeStart(const String & filename) nogil except + 
            # wrap-doc:
            #   Prepares a FASTA file given by 'filename' for streamed writing using writeNext()
            #   -----
            #   :raises:
            #       Exception:UnableToCreateFile is thrown if the process is not able to write to the file (disk full?)
        void writeNext(const FASTAEntry & protein) nogil except +
            # wrap-doc:
            #   Stores the data given by `protein`. Call writeStart() once before calling writeNext()
            #   -----
            #   Call writeEnd() when done to close the file!
            #   -----
            #   :raises:
            #       Exception:UnableToCreateFile is thrown if the process is not able to write to the file (disk full?)
        void writeEnd() nogil except + # wrap-doc:Closes the file (flush). Called implicitly when FASTAFile object does out of scope

cdef extern from "<OpenMS/FORMAT/FASTAFile.h>" namespace "OpenMS::FASTAFile":

    cdef cppclass FASTAEntry:
        FASTAEntry() nogil except +
        FASTAEntry(FASTAEntry) nogil except +

        String identifier
        String description
        String sequence

        bool headerMatches(const FASTAEntry & rhs) nogil except +
        bool sequenceMatches(const FASTAEntry & rhs) nogil except +
