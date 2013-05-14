
from MSExperiment  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/MzDataFile.h>" namespace "OpenMS":

    cdef cppclass MzDataFile(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        MzDataFile() nogil except +

        void load(String, MSExperiment[Peak1D, ChromatogramPeak] &) nogil except+
        void store(String, MSExperiment[Peak1D, ChromatogramPeak] &) nogil except+

        PeakFileOptions getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +
