
from MSExperiment  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from StringList cimport *
from ProgressLogger cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/MzDataFile.h>" namespace "OpenMS":

    cdef cppclass MzDataFile(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        MzDataFile() nogil except +

        void load(const String&, MSExperiment &) nogil except+
        void store(const String&, MSExperiment &) nogil except+

        PeakFileOptions getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +

        bool isSemanticallyValid(const String& filename, StringList & errors, StringList & warnings) nogil except +

