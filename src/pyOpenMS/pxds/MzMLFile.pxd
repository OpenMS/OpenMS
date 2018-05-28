from MSExperiment  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from PeakFileOptions cimport *
from IMSDataConsumer cimport *

cdef extern from "<OpenMS/FORMAT/MzMLFile.h>" namespace "OpenMS":

    cdef cppclass MzMLFile(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        MzMLFile() nogil except +

        void load(const String&, MSExperiment &) nogil except+
        void store(const String&, MSExperiment &) nogil except+

        void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *) nogil except + # wrap-ignore

        PeakFileOptions getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +

        bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings) nogil except +

        # NAMESPACE # bool isValid(const String & filename, std::ostream & os)

        # void loadSize(const String & filename, Size & scount, Size & ccount) 
