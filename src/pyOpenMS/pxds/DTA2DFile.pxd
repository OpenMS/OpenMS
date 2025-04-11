from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from String cimport *
from ProgressLogger cimport *

from MSRun cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/DTA2DFile.h>" namespace "OpenMS":

    cdef cppclass DTA2DFile(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        DTA2DFile() except + nogil 
        DTA2DFile(DTA2DFile &) except + nogil  # compiler

        void storeTIC(String filename, MSRun & peakmap) except + nogil 
        void store(String filename, MSRun & peakmap) except + nogil 
        void load(String filename, MSRun & peakmap) except + nogil 
        PeakFileOptions  getOptions() except + nogil 

