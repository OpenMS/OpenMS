from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from String cimport *
from ProgressLogger cimport *

from MSSpectrum cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FORMAT/DTAFile.h>" namespace "OpenMS":

    cdef cppclass DTAFile:

        DTAFile() except + nogil 
        DTAFile(DTAFile &) except + nogil  # compiler

        void load(String filename, MSSpectrum & spectrum) except + nogil 
        void store(String filename, MSSpectrum & spectrum) except + nogil 
