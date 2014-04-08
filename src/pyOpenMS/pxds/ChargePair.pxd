from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from Compomer cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/ChargePair.h>" namespace "OpenMS":

    cdef cppclass ChargePair:
  
        ChargePair() nogil except +
        ChargePair(ChargePair) nogil except + 
  
        ChargePair(Size index0,
               Size index1,
               Int charge0,
               Int charge1,
               Compomer compomer,
               double mass_diff,
                   bool active) nogil except +

        Int getCharge(UInt pairID) nogil except +
    
        void setCharge(UInt pairID, Int e) nogil except +
    
        Size getElementIndex(UInt pairID) nogil except +
    
        void setElementIndex(UInt pairID, Size e) nogil except +
    
        Compomer getCompomer() nogil except +
    
        void setCompomer( Compomer & compomer) nogil except +
    
        double getMassDiff() nogil except +
    
        void setMassDiff(double mass_diff) nogil except +
    
        double getEdgeScore() nogil except +
    
        void setEdgeScore(double score) nogil except +
    
        bool isActive() nogil except +
    
        void setActive( bool active) nogil except +
