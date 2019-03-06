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

        void load(const String& filename, MSExperiment &) nogil except+
        void store(const String& filename, MSExperiment &) nogil except+

        # COMMENT: store/load XML structure to/from a string
        void storeBuffer(String & output, MSExperiment exp) nogil except +
        void loadBuffer(const String& input, MSExperiment & exp) nogil except +

        void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *) nogil except + # wrap-ignore
        void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *,
                       bool skip_full_count, bool skip_first_pass) nogil except + # wrap-ignore

        void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, MSExperiment& e) nogil except + # wrap-ignore
        void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, MSExperiment& e,
                       bool skip_full_count, bool skip_first_pass) nogil except + # wrap-ignore

        PeakFileOptions getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +

        bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings) nogil except +

