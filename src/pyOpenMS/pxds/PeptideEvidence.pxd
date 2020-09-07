from Types cimport *
from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/METADATA/PeptideEvidence.h>" namespace "OpenMS":
    cdef cppclass PeptideEvidence :

        PeptideEvidence() nogil except +
        PeptideEvidence(PeptideEvidence) nogil except +

        # const members
        ## int UNKNOWN_POSITION
        ## int N_TERMINAL_POSITION
        ## char UNKNOWN_AA
        ## char N_TERMINAL_AA
        ## char C_TERMINAL_AA

        void setStart(Int start) nogil except +
        Int getStart() nogil except +
        void setEnd(Int end) nogil except +
        Int getEnd() nogil except +
        void setAABefore(char rhs) nogil except +
        char getAABefore() nogil except +
        void setAAAfter(char rhs) nogil except +
        char getAAAfter() nogil except +
        void setProteinAccession(String s) nogil except +
        String getProteinAccession() nogil except +

        bool hasValidLimits() nogil except +

        bool operator==(PeptideEvidence & rhs) nogil except +
        bool operator!=(PeptideEvidence & rhs) nogil except +

