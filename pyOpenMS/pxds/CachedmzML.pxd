from MSExperiment  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h>" namespace "OpenMS":

    cdef cppclass CachedmzML(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger


        CachedmzML() nogil except +
        CachedmzML(CachedmzML) nogil except +

        void writeMemdump(MSExperiment[Peak1D,ChromatogramPeak] exp, String out) nogil except +
        void writeMetadata(MSExperiment[Peak1D,ChromatogramPeak] exp, String out_meta) nogil except +

        void readMemdump(MSExperiment[Peak1D,ChromatogramPeak] exp, String filename) nogil except +

