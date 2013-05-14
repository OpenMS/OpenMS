from Types cimport *
from String cimport *

cdef extern from "<OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>" namespace "OpenMS":

    cdef cppclass ReactionMonitoringTransition:

        ReactionMonitoringTransition()                  nogil except +
        ReactionMonitoringTransition(ReactionMonitoringTransition)   nogil except + #wrap-ignore
        String getName() 
        String getNativeID() 
        String getPeptideRef() 
        void setName(String & name) 
        void setNativeID(String & name) 
        void setPeptideRef(String & peptide_ref) 

        DoubleReal getProductMZ() except +
        void setProductMZ(DoubleReal)

        DoubleReal getPrecursorMZ() except +
        void setPrecursorMZ(DoubleReal)

