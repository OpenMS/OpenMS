from SourceFile cimport *
from DateTime cimport *
from FileTypes cimport *
from String cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/METADATA/DocumentIdentifier.h>" namespace "OpenMS":

    cdef cppclass DocumentIdentifier:

        DocumentIdentifier() nogil except +
        DocumentIdentifier(DocumentIdentifier) nogil except + # wrap-ignore

        # set document identifier (e.g. an LSID)
        void setIdentifier(String id) nogil except + # wrap-doc:Set document identifier (e.g. an LSID)

        # retrieve document identifier (e.g. an LSID)
        String getIdentifier() nogil except + # wrap-doc:Retrieve document identifier (e.g. an LSID)

        # exchange content with @p from -> gets overwritten later
        # this does not work since some derived classes overwrite this (e.g. MSExperiment)
        ## void swap(DocumentIdentifier from_) nogil except +

        # set the file_type according to the type of the file loaded from (see FileHandler::Type) preferably done whilst loading
        void setLoadedFileType(String file_name) nogil except + # wrap-doc:Set the file_type according to the type of the file loaded from, preferably done whilst loading

        # get the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded from
        int getLoadedFileType() nogil except + # wrap-doc:Get the file_type (e.g. featureXML, consensusXML, mzData, mzXML, mzML, ...) of the file loaded from

        # set the file_name_ according to absolute path of the file loaded from preferably done whilst loading
        void setLoadedFilePath(String file_name) nogil except + # wrap-doc:Set the file_name_ according to absolute path of the file loaded from preferably done whilst loading
        # get the file_name_ which is the absolute path to the file loaded from
        String getLoadedFilePath() nogil except + # wrap-doc:Get the file_name which is the absolute path to the file loaded

        # Errors in C++
        # void swap(DocumentIdentifier & from_) nogil except +
