from MSExperiment  cimport *
from FeatureMap cimport *
from Feature cimport *
from libcpp.string cimport string as libcpp_string
from FileTypes cimport *
from Types cimport *

cdef extern from "<OpenMS/FORMAT/FileHandler.h>" namespace "OpenMS":

    cdef cppclass FileHandler:  # wrap=True
        FileHandler() nogil except +
        FileHandler(FileHandler) nogil except +

        void loadExperiment(libcpp_string, MSExperiment[Peak1D, ChromatogramPeak] &) nogil except+
        void storeExperiment(libcpp_string, MSExperiment[Peak1D, ChromatogramPeak]) nogil except+
        void loadFeatures(libcpp_string, FeatureMap[Feature] &) nogil except +



cdef extern from "<OpenMS/FORMAT/FileHandler.h>" namespace "OpenMS::FileHandler":

    int  getType(String filename) nogil except + # wrap-attach:FileHandler
