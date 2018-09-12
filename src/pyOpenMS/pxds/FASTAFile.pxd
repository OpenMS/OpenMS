from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool

from String cimport *
from DefaultParamHandler cimport *

from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/FORMAT/FASTAFile.h>" namespace "OpenMS":

    cdef cppclass FASTAFile:

        FASTAFile() nogil except +
        FASTAFile(FASTAFile) nogil except + #wrap-ignore

        void load(const String& filename, libcpp_vector[FASTAEntry] & data) nogil except +
        void store(const String& filename, libcpp_vector[FASTAEntry] & data) nogil except +

        void readStart(const String & filename) nogil except +
        bool readNext(FASTAEntry & protein) nogil except +
        # NAMESPACE # std::streampos position() nogil except +
        bool atEnd() nogil except +
        # NAMESPACE # bool setPosition(const std::streampos & pos) nogil except +
        void writeStart(const String & filename) nogil except +
        void writeNext(const FASTAEntry & protein) nogil except +
        void writeEnd() nogil except +

cdef extern from "<OpenMS/FORMAT/FASTAFile.h>" namespace "OpenMS::FASTAFile":

    cdef cppclass FASTAEntry:
        FASTAEntry() nogil except +
        FASTAEntry(FASTAEntry) nogil except +

        String identifier
        String description
        String sequence

        bool headerMatches(const FASTAEntry & rhs) nogil except +
        bool sequenceMatches(const FASTAEntry & rhs) nogil except +
