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

        void load(String, MSExperiment[Peak1D, ChromatogramPeak] &) nogil except+
        void store(String, MSExperiment[Peak1D, ChromatogramPeak] &) nogil except+

        void transform(String, IMSDataConsumer[Peak1D, ChromatogramPeak] *) nogil except + # wrap-ignore

        PeakFileOptions getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +

        bool isSemanticallyValid(String & filename, StringList & errors, StringList & warnings) nogil except +

        # NAMESPACE # bool isValid(String & filename, std::ostream & os)

        # void loadSize(String & filename, Size & scount, Size & ccount) 
