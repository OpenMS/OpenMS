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
        #
        # wrap-doc:
        #   File adapter for MzML files
        #   -----
        #   Provides methods to load and store MzML files.
        #   PeakFileOptions allow to load a reduced subset of the data into an MSExperiment.
        #   -----
        #   See help(MSExperiment) how data is stored after loading.
        #   See help(PeakFileOptions) for available options.
        #   -----
        #   Usage:
        #     exp = MSExperiment()
        #     MzMLFile().load("test.mzML", exp)
        #     spec = []
        #     for s in exp.getSpectra():
        #       if s.getMSLevel() != 1:
        #         spec.append(s)
        #     exp.setSpectra(spec)
        #     MzMLFile().store("filtered.mzML", exp)
        #   -----        

        MzMLFile() nogil except +

        void load(const String& filename, MSExperiment &) nogil except+ #wrap-doc: Loads from an MzML file. Spectra and chromatograms are sorted by default (this can be disabled using PeakFileOptions).
        void store(const String& filename, MSExperiment &) nogil except+ #wrap-doc: Stores a map in an MzML file.

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
        void setOptions(PeakFileOptions) nogil except + #wrap-doc: set PeakFileOptions to perform filtering during loading. E.g., to load only MS1 spectra or meta data only.

        bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings) nogil except +

