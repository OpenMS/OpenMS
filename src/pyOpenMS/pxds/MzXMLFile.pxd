from MSExperiment  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from PeakFileOptions cimport *
from IMSDataConsumer cimport *

cdef extern from "<OpenMS/FORMAT/MzXMLFile.h>" namespace "OpenMS":

    cdef cppclass MzXMLFile:
        # wrap-doc:
        #  File adapter for MzXML files
        #  
        #  Provides methods to load and store MzXML files.
        #  MzXML is an older format; for new projects consider using MzML instead.
        #  
        #  Usage:
        #  
        #  .. code-block:: python
        #  
        #    exp = MSExperiment()
        #    MzXMLFile().load("test.mzXML", exp)

        MzXMLFile() except + nogil
        MzXMLFile(MzXMLFile &) except + nogil  #compiler

        void load(String filename, MSExperiment & exp) except + nogil
            # wrap-doc:
                #  Loads a MSExperiment from a MzXML file
                #  
                #  
                #  :param exp: MSExperiment

        void store(String filename, MSExperiment & exp) except + nogil
            # wrap-doc:
                #  Stores a MSExperiment in a MzXML file
                #  
                #  
                #  :param exp: MSExperiment

        void transform(String, IMSDataConsumer[Peak1D, ChromatogramPeak] *) except + nogil  # wrap-ignore

        PeakFileOptions getOptions() except + nogil  # wrap-doc:Returns the options for loading/storing
        void setOptions(PeakFileOptions) except + nogil  # wrap-doc:Sets options for loading/storing
