from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from ProteinIdentification cimport *
from PeptideIdentification cimport *

from String cimport *
from ProgressLogger cimport *
from Identification cimport *

cdef extern from "<OpenMS/FORMAT/MzIdentMLFile.h>" namespace "OpenMS":

    cdef cppclass MzIdentMLFile(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        MzIdentMLFile() nogil except +

        void load(String filename, libcpp_vector[ProteinIdentification] & poid, libcpp_vector[PeptideIdentification] & peid) nogil except +
        void store(String filename, libcpp_vector[ProteinIdentification] & poid, libcpp_vector[PeptideIdentification] & peid) nogil except +
        bool isSemanticallyValid(String filename, StringList errors, StringList warnings) nogil except +

