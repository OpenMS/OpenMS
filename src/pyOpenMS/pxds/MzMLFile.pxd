from MSRun  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from PeakFileOptions cimport *
from IMSDataConsumer cimport *

cdef extern from "<OpenMS/FORMAT/MzMLFile.h>" namespace "OpenMS":

    cdef cppclass MzMLFile(ProgressLogger):
        # wrap-inherits:
        #  ProgressLogger
        #
        # wrap-doc:
        #  File adapter for MzML files
        #  
        #  Provides methods to load and store MzML files.
        #  PeakFileOptions allow to load a reduced subset of the data into an MSRun.
        #  
        #  See help(MSRun) how data is stored after loading.
        #  See help(PeakFileOptions) for available options.
        #  
        #  Usage:
        #  
        #  .. code-block:: python
        #  
        #    exp = MSRun()
        #    MzMLFile().load("test.mzML", exp)
        #    spec = []
        #    for s in exp.getSpectra():
        #      if s.getMSLevel() != 1:
        #        spec.append(s)
        #    exp.setSpectra(spec)
        #    MzMLFile().store("filtered.mzML", exp)
        #          

        MzMLFile() except + nogil 
        MzMLFile(MzMLFile &) except + nogil 

        void load(const String& filename, MSRun &) except + nogil # wrap-doc:Loads from an MzML file. Spectra and chromatograms are sorted by default (this can be disabled using PeakFileOptions)
        void store(const String& filename, MSRun &) except + nogil # wrap-doc:Stores a MSRun in an MzML file

        # COMMENT: store/load XML structure to/from a string
        void storeBuffer(String & output, MSRun exp) except + nogil 
            # wrap-doc:
                #  Stores a map in an output string
                #  
                #  
                #  :param output: An empty string to store the result
                #  :param exp: Has to be an MSRun

        void loadBuffer(const String& input, MSRun & exp) except + nogil 
            # wrap-doc:
                #  Loads a map from a MzML file stored in a buffer (in memory)
                #  
                #  
                #  :param buffer: The buffer with the data (i.e. string with content of an mzML file)
                #  :param exp: Is an MSRun
                #  :raises:
                #    Exception: ParseError is thrown if an error occurs during parsing

        void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *) except + nogil  # wrap-ignore
        void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *,
                       bool skip_full_count, bool skip_first_pass) except + nogil  # wrap-ignore

        void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, MSRun& e) except + nogil  # wrap-ignore
        void transform(const String&, IMSDataConsumer[Peak1D, ChromatogramPeak] *, MSRun& e,
                       bool skip_full_count, bool skip_first_pass) except + nogil  # wrap-ignore

        PeakFileOptions getOptions() except + nogil 
        void setOptions(PeakFileOptions) except + nogil  # wrap-doc:Set PeakFileOptions to perform filtering during loading. E.g., to load only MS1 spectra or meta data only

        bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings) except + nogil 
