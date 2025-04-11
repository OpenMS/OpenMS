from MSRun  cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from String cimport *
from ProgressLogger cimport *
from PeakFileOptions cimport *
from IMSDataConsumer cimport *

cdef extern from "<OpenMS/FORMAT/MzXMLFile.h>" namespace "OpenMS":

    cdef cppclass MzXMLFile(ProgressLogger):
        # wrap-inherits:
        #  ProgressLogger

        MzXMLFile() except + nogil 
        MzXMLFile(MzXMLFile &) except + nogil  #compiler

        void load(String filename, MSRun & exp) except + nogil  
            # wrap-doc:
                #  Loads a MSRun from a MzXML file
                #  
                #  
                #  :param exp: MSRun

        void store(String filename, MSRun & exp) except + nogil 
            # wrap-doc:
                #  Stores a MSRun in a MzXML file
                #  
                #  
                #  :param exp: MSRun

        void transform(String, IMSDataConsumer[Peak1D, ChromatogramPeak] *) except + nogil  # wrap-ignore

        PeakFileOptions getOptions() except + nogil  # wrap-doc:Returns the options for loading/storing
        void setOptions(PeakFileOptions) except + nogil  # wrap-doc:Sets options for loading/storing
