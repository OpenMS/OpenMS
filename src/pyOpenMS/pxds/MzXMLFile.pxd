from MSExperiment  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from PeakFileOptions cimport *
from IMSDataConsumer cimport *

cdef extern from "<OpenMS/FORMAT/MzXMLFile.h>" namespace "OpenMS":

    cdef cppclass MzXMLFile(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        MzXMLFile() nogil except +
        MzXMLFile(MzXMLFile &) nogil except + #compiler

        void load(String filename, MSExperiment & exp) nogil except + 
            # wrap-doc:
                #   Loads a MSExperiment from a MzXML file
                #   -----
                #   :param exp: MSExperiment

        void store(String filename, MSExperiment & exp) nogil except +
            # wrap-doc:
                #   Stores a MSExperiment in a MzXML file
                #   -----
                #   :param exp: MSExperiment

        void transform(String, IMSDataConsumer[Peak1D, ChromatogramPeak] *) nogil except + # wrap-ignore

        PeakFileOptions getOptions() nogil except + # wrap-doc:Returns the options for loading/storing
        void setOptions(PeakFileOptions) nogil except + # wrap-doc:Sets options for loading/storing
