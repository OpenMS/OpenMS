from MSExperiment  cimport *
from FeatureMap cimport *
from Feature cimport *
from String cimport *
from libcpp.string cimport string as libcpp_string
from FileTypes cimport *
from Types cimport *
from PeakFileOptions cimport *

cdef extern from "<OpenMS/FORMAT/FileHandler.h>" namespace "OpenMS":
        # wrap-doc:
        #   Facilitates file handling by file type recognition.
        #   This class provides file type recognition from the file name and
        #   for some types from the file content.
        #   It offers a common interface to load MSExperiment data
        #   and allows querying for supported file types.
        #   -----
        #   Usage:
        #     MSExperiment exp;
        #     FileHandler().loadExperiment("test.mzXML", exp)
        #     FileHandler().loadExperiment("test.mzML", exp)
        #   -----

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

    # Returns the file name without the extension
    String stripExtension(String file) nogil except + # wrap-attach:FileHandler
    # Removes the current extension (if any) and adds a new one
    String swapExtension(String filename, FileType new_type) nogil except + # wrap-attach:FileHandler 
