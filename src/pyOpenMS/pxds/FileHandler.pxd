from MSExperiment  cimport *
from FeatureMap cimport *
from Feature cimport *
from String cimport *
from libcpp.string cimport string as libcpp_string
from FileTypes cimport *
from Types cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/FileHandler.h>" namespace "OpenMS":

    cdef cppclass FileHandler:  # wrap=True
        FileHandler() nogil except +
        FileHandler(FileHandler) nogil except +

        bool loadExperiment(String, MSExperiment &) nogil except+
        void storeExperiment(String, MSExperiment) nogil except+
        bool loadFeatures(String, FeatureMap &) nogil except +

        PeakFileOptions  getOptions() nogil except +
        void setOptions(PeakFileOptions) nogil except +

#
# wrap static method:
#
cdef extern from "<OpenMS/FORMAT/FileHandler.h>" namespace "OpenMS::FileHandler":

    int getType(const String& filename) nogil except + # wrap-attach:FileHandler
    FileType getTypeByFileName(const String & filename) nogil except + # wrap-attach:FileHandler 
    FileType getTypeByContent(const String & filename) nogil except + # wrap-attach:FileHandler 
    String computeFileHash(const String & filename) nogil except + # wrap-attach:FileHandler 
    bool isSupported(FileType type_) nogil except + # wrap-attach:FileHandler 
    bool hasValidExtension(const String & filename, FileType type_) nogil except + # wrap-attach:FileHandler 

