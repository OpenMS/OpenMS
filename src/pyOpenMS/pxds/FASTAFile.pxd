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


cdef extern from "<OpenMS/FORMAT/FASTAFile.h>" namespace "OpenMS::FASTAFile":

    cdef cppclass FASTAEntry:
        FASTAEntry() nogil except +
        FASTAEntry(FASTAEntry) nogil except +

        String identifier
        String description
        String sequence

