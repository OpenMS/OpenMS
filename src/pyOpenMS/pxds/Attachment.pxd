from libcpp cimport bool
from ProgressLogger cimport *
from UniqueIdGenerator cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/QcMLFile.h>" namespace "OpenMS::QcMLFile":
    
    cdef cppclass Attachment "OpenMS::QcMLFile::Attachment":
        Attachment() nogil except +
        Attachment(Attachment) nogil except +
        String name
        String value
        String cvRef
        String cvAcc
        String unitRef
        String unitAcc
        String binary
        String qualityRef
        libcpp_vector[ String ] colTypes
        libcpp_vector[ libcpp_vector[ String ] ] tableRows # wrap-ignore
        bool operator==(Attachment &rhs) nogil except +
        bool operator<(Attachment &rhs) nogil except +
        bool operator>(Attachment &rhs) nogil except +
        String toXMLString(UInt indentation_level) nogil except +
        String toCSVString(String separator) nogil except +

