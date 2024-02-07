from SourceFile cimport *
from DateTime cimport *
from FileTypes cimport *
from String cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/DocumentIdentifier.h>" namespace "OpenMS":

    cdef cppclass DocumentIdentifier:

        DocumentIdentifier() except + nogil 
        DocumentIdentifier(DocumentIdentifier &) except + nogil 

        void setIdentifier(String id) except + nogil  # wrap-doc:Sets document identifier (e.g. an LSID)

        String getIdentifier() except + nogil  # wrap-doc:Retrieve document identifier (e.g. an LSID)

        # exchange content with @p from -> gets overwritten later
        # this does not work since some derived classes overwrite this (e.g. MSExperiment)
        ## void swap(DocumentIdentifier from_) except + nogil 

        void setLoadedFileType(String file_name) except + nogil  # wrap-doc:Sets the file_type according to the type of the file loaded from, preferably done whilst loading

        int getLoadedFileType() except + nogil  # wrap-doc:Returns the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded

        void setLoadedFilePath(String file_name) except + nogil  # wrap-doc:Sets the file_name according to absolute path of the file loaded, preferably done whilst loading
        String getLoadedFilePath() except + nogil  # wrap-doc:Returns the file_name which is the absolute path to the file loaded

        # Errors in C++
        # void swap(DocumentIdentifier & from_) except + nogil 
