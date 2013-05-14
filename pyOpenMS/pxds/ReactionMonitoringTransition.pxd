from Types cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS":

    cdef cppclass ReactionMonitoringTransition:

        ReactionMonitoringTransition()                  nogil except +
        ReactionMonitoringTransition(ReactionMonitoringTransition)   nogil except + #wrap-ignore
        String getName()                           nogil except +
        String getNativeID()                       nogil except +
        String getPeptideRef()                     nogil except +
        void setName(String name)                  nogil except +
        void setNativeID(String name)              nogil except +
        void setPeptideRef(String peptide_ref)     nogil except +

        DoubleReal getProductMZ()         nogil except +
        void setProductMZ(DoubleReal)     nogil except +

        DoubleReal getPrecursorMZ()       nogil except +
        void setPrecursorMZ(DoubleReal)   nogil except +
      
        DecoyTransitionType getDecoyTransitionType() nogil except +


cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS::ReactionMonitoringTransition":

    cdef enum DecoyTransitionType:

        UNKNOWN, TARGET, DECOY

