from MSExperiment  cimport *
from FeatureMap cimport *
from Feature cimport *
from libcpp.string cimport string as libcpp_string
from FileTypes cimport *
from Types cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/FileHandler.h>" namespace "OpenMS":

    cdef cppclass FileHandler:  # wrap=True
        FileHandler() nogil except +
        FileHandler(FileHandler) nogil except +

        bool loadExperiment(libcpp_string, MSExperiment &) nogil except+
        void storeExperiment(libcpp_string, MSExperiment) nogil except+
        bool loadFeatures(libcpp_string, FeatureMap &) nogil except +

        PeakFileOptions  getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +

#
# wrap static method:
#
cdef extern from "<OpenMS/FORMAT/FileHandler.h>" namespace "OpenMS::FileHandler":

    int  getType(String filename) nogil except + # wrap-attach:FileHandler
    FileType getTypeByFileName(String & filename) nogil except + # wrap-attach:FileHandler 
    FileType getTypeByContent(String & filename) nogil except + # wrap-attach:FileHandler 
    String computeFileHash(String & filename) nogil except + # wrap-attach:FileHandler 
    bool isSupported(FileType type_) nogil except + # wrap-attach:FileHandler 
