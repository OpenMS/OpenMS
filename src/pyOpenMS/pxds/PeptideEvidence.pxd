from Types cimport *
from libcpp cimport bool
from Types cimport *
from String cimport *
from MetaInfoInterface cimport *
from AASequence cimport *

cdef extern from "<OpenMS/METADATA/PeptideEvidence.h>" namespace "OpenMS":
    
    cdef cppclass PeptideEvidence(MetaInfoInterface) :
        # wrap-inherits:
        #  MetaInfoInterface
        PeptideEvidence() nogil except +
        PeptideEvidence(PeptideEvidence) nogil except +
        void setStart(Int start) nogil except +
        Int getStart() nogil except +
        void setEnd(Int end) nogil except +
        Int getEnd() nogil except +
        void setAABefore(char rhs) nogil except +
        char getAABefore() nogil except +
        void setAAAfter(char rhs) nogil except +
        char getAAAfter() nogil except +
        bool operator==(PeptideEvidence & rhs) nogil except +
        bool operator!=(PeptideEvidence & rhs) nogil except +

