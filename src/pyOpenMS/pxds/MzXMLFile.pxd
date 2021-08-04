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

        void load(String, MSExperiment &) nogil except + 
            # wrap-doc:
                #   Loads a map from a MzXML file
                #   -----
                #   :param map: Has to be a MSExperiment or have the same interface

        void store(String, MSExperiment &) nogil except +
            # wrap-doc:
                #   Stores a map in a MzXML file
                #   -----
                #   :param map: Has to be a MSExperiment or have the same interface

        void transform(String, IMSDataConsumer[Peak1D, ChromatogramPeak] *) nogil except + # wrap-ignore

        PeakFileOptions getOptions() nogil except + # wrap-doc:Returns the options for loading/storing
        void setOptions(PeakFileOptions) nogil except + # wrap-doc:Set options for loading/storing
