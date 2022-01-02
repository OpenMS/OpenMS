from SourceFile cimport *
from DateTime cimport *
from FileTypes cimport *
from String cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/DocumentIdentifier.h>" namespace "OpenMS":

    cdef cppclass DocumentIdentifier:

        DocumentIdentifier() nogil except +
        DocumentIdentifier(DocumentIdentifier &) nogil except +

        void setIdentifier(String id) nogil except + # wrap-doc:Sets document identifier (e.g. an LSID)

        String getIdentifier() nogil except + # wrap-doc:Retrieve document identifier (e.g. an LSID)

        # exchange content with @p from -> gets overwritten later
        # this does not work since some derived classes overwrite this (e.g. MSExperiment)
        ## void swap(DocumentIdentifier from_) nogil except +

        void setLoadedFileType(String file_name) nogil except + # wrap-doc:Sets the file_type according to the type of the file loaded from, preferably done whilst loading

        int getLoadedFileType() nogil except + # wrap-doc:Returns the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded

        void setLoadedFilePath(String file_name) nogil except + # wrap-doc:Sets the file_name according to absolute path of the file loaded, preferably done whilst loading
        String getLoadedFilePath() nogil except + # wrap-doc:Returns the file_name which is the absolute path to the file loaded

        # Errors in C++
        # void swap(DocumentIdentifier & from_) nogil except +
