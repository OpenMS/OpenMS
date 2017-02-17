from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool

from String cimport *
from DefaultParamHandler cimport *

from PeptideIdentification cimport *
from ProteinIdentification cimport *

cdef extern from "<OpenMS/FORMAT/FASTQFile.h>" namespace "OpenMS":

    cdef cppclass FASTQFile:

        FASTQFile() nogil except +
        FASTQFile(FASTQFile) nogil except + #wrap-ignore

        void load(String& filename, libcpp_vector[FASTQEntry] & data) nogil except +
        void store(String& filename, libcpp_vector[FASTQEntry] & data) nogil except +


cdef extern from "<OpenMS/FORMAT/FASTQFile.h>" namespace "OpenMS::FASTQFile":

    cdef cppclass FASTQEntry:
        FASTQEntry() nogil except +
        FASTQEntry(FASTQEntry) nogil except +

        String identifier
        String description
        String sequence
		String quality
		