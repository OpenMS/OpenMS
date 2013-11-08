from Types  cimport *

from PeakFileOptions cimport *
from MzMLFile cimport *
from OnDiscMSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/IndexedMzMLFileLoader.h>" namespace "OpenMS":

    cdef cppclass IndexedMzMLFileLoader:

        IndexedMzMLFileLoader() nogil except +

        void load(String, OnDiscMSExperiment[Peak1D, ChromatogramPeak] &) nogil except+
        void store(String, OnDiscMSExperiment[Peak1D, ChromatogramPeak] &) nogil except+
        void store(String, MSExperiment[Peak1D, ChromatogramPeak] &) nogil except+

        PeakFileOptions getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +

