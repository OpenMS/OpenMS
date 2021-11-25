from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from Compomer cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/ChargePair.h>" namespace "OpenMS":

    cdef cppclass ChargePair:
  
        ChargePair() nogil except +
        ChargePair(ChargePair &) nogil except +
  
        ChargePair(Size index0,
               Size index1,
               Int charge0,
               Int charge1,
               Compomer compomer,
               double mass_diff,
               bool active) nogil except +

        Int getCharge(UInt pairID) nogil except + # wrap-doc:Returns the charge (for element 0 or 1)
        void setCharge(UInt pairID, Int e) nogil except + # wrap-doc:Sets the charge (for element 0 or 1)
    
        Size getElementIndex(UInt pairID) nogil except + # wrap-doc:Returns the element index (for element 0 or 1)
        void setElementIndex(UInt pairID, Size e) nogil except + # wrap-doc:Sets the element index (for element 0 or 1)
    
        Compomer getCompomer() nogil except + # wrap-doc:Returns the Id of the compomer that explains the mass difference
        void setCompomer( Compomer & compomer) nogil except + # wrap-doc:Sets the compomer id
    
        double getMassDiff() nogil except + # wrap-doc:Returns the mass difference
        void setMassDiff(double mass_diff) nogil except + # wrap-doc:Sets the mass difference
    
        double getEdgeScore() nogil except + # wrap-doc:Returns the ILP edge score
        void setEdgeScore(double score) nogil except + # wrap-doc:Sets the ILP edge score
    
        bool isActive() nogil except + # wrap-doc:Is this pair realized?
        void setActive( bool active) nogil except +
