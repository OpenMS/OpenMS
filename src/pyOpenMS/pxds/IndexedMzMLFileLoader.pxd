from Types cimport *

from PeakFileOptions cimport *
from MzMLFile cimport *
from OnDiscMSExperiment cimport *

cdef extern from "<OpenMS/FORMAT/IndexedMzMLFileLoader.h>" namespace "OpenMS":

    cdef cppclass IndexedMzMLFileLoader:

        IndexedMzMLFileLoader() nogil except +
 
        bool load(String, OnDiscMSExperiment &) nogil except+
        void store(String, OnDiscMSExperiment &) nogil except+
        void store(String, MSExperiment &) nogil except+

        PeakFileOptions getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +

