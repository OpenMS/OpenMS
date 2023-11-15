from libcpp cimport bool
from ProgressLogger cimport *
from UniqueIdGenerator cimport *
from String cimport *

cdef extern from "<OpenMS/FORMAT/QcMLFile.h>" namespace "OpenMS::QcMLFile":
    
    cdef cppclass Attachment "OpenMS::QcMLFile::Attachment":
        Attachment() except + nogil 
        Attachment(Attachment &) except + nogil 
        String name
        String id
        String value
        String cvRef
        String cvAcc
        String unitRef
        String unitAcc
        String binary
        String qualityRef
        libcpp_vector[ String ] colTypes
        libcpp_vector[ libcpp_vector[ String ] ] tableRows # wrap-ignore
        bool operator==(Attachment &rhs) except + nogil 
        bool operator<(Attachment &rhs) except + nogil 
        bool operator>(Attachment &rhs) except + nogil 
        String toXMLString(UInt indentation_level) except + nogil 
        String toCSVString(String separator) except + nogil 

