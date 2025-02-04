from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from Types cimport *
from Compomer cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/ChargePair.h>" namespace "OpenMS":

    cdef cppclass ChargePair:
  
        ChargePair() except + nogil 
        ChargePair(ChargePair &) except + nogil 
  
        ChargePair(Size index0,
               Size index1,
               Int charge0,
               Int charge1,
               Compomer compomer,
               double mass_diff,
               bool active) except + nogil 

        Int getCharge(UInt pairID) except + nogil  # wrap-doc:Returns the charge (for element 0 or 1)
        void setCharge(UInt pairID, Int e) except + nogil  # wrap-doc:Sets the charge (for element 0 or 1)
    
        Size getElementIndex(UInt pairID) except + nogil  # wrap-doc:Returns the element index (for element 0 or 1)
        void setElementIndex(UInt pairID, Size e) except + nogil  # wrap-doc:Sets the element index (for element 0 or 1)
    
        Compomer getCompomer() except + nogil  # wrap-doc:Returns the Id of the compomer that explains the mass difference
        void setCompomer( Compomer & compomer) except + nogil  # wrap-doc:Sets the compomer id
    
        double getMassDiff() except + nogil  # wrap-doc:Returns the mass difference
        void setMassDiff(double mass_diff) except + nogil  # wrap-doc:Sets the mass difference
    
        double getEdgeScore() except + nogil  # wrap-doc:Returns the ILP edge score
        void setEdgeScore(double score) except + nogil  # wrap-doc:Sets the ILP edge score
    
        bool isActive() except + nogil  # wrap-doc:Is this pair realized?
        void setActive( bool active) except + nogil 
