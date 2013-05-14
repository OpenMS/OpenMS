from Types cimport *
from libcpp cimport bool
from Types cimport *
from String cimport *
from MetaInfoInterface cimport *
# TODO add in OpenMS the import
from AASequence cimport *

cdef extern from "<OpenMS/METADATA/PeptideEvidence.h>" namespace "OpenMS":
    
    cdef cppclass PeptideEvidence(MetaInfoInterface) :
        # wrap-inherits:
        #  MetaInfoInterface
        PeptideEvidence() nogil except +
        PeptideEvidence(PeptideEvidence) nogil except +
        # TODO this constructor does not exist, remove it from OpenMS
        # PeptideEvidence(DoubleReal score, UInt rank, Int charge, AASequence & sequence) nogil except +
        String  getDBSequenceRef() nogil except +
        void setDBSequenceRef(String & rhs) nogil except +
        String  getTranslationTableRef() nogil except +
        void setTranslationTableRef(String & rhs) nogil except +
        void setStart(Int start) nogil except +
        Int getStart() nogil except +
        void setEnd(Int end) nogil except +
        Int getEnd() nogil except +
        void setPre(char rhs) nogil except +
        char getPre() nogil except +
        void setPost(char rhs) nogil except +
        char getPost() nogil except +
        void setId(String & id_) nogil except +
        String  getId() nogil except +
        void setName(String & name) nogil except +
        String  getName() nogil except +
        void setMissedCleavages(Int rhs) nogil except +
        Int getMissedCleavages() nogil except +
        void setIsDecoy(bool is_decoy) nogil except +
        bool getIsDecoy() nogil except +
        void setFrame(Int frame) nogil except +
        Int getFrame() nogil except +
        bool operator==(PeptideEvidence & rhs) nogil except +
        bool operator!=(PeptideEvidence & rhs) nogil except +

