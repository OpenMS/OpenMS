from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool

from String cimport *
from ProgressLogger cimport *

from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/DTA2DFile.h>" namespace "OpenMS":

    cdef cppclass DTA2DFile(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        DTA2DFile() nogil except +
        DTA2DFile(DTA2DFile &) nogil except + # compiler

        void storeTIC(String filename, MSExperiment & peakmap) nogil except +
        void store(String filename, MSExperiment & peakmap) nogil except +
        void load(String filename, MSExperiment & peakmap) nogil except +
        PeakFileOptions  getOptions() nogil except +

