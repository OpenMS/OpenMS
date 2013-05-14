from libcpp.vector cimport vector as libcpp_vector
from String cimport *
from Peak1D cimport *
from RichPeak1D cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/MSPFile.h>" namespace "OpenMS":

    cdef cppclass MSPFile:

        MSPFile() nogil except +

        void store(String filename, MSExperiment[RichPeak1D, ChromatogramPeak] & exp)
        void load(String filename, libcpp_vector[PeptideIdentification] & ids, MSExperiment[RichPeak1D, ChromatogramPeak] & exp)
