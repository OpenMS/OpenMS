from libcpp.vector cimport vector as libcpp_vector
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/IsotopeCluster.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeCluster "OpenMS::IsotopeCluster":

        IsotopeCluster() nogil except + # wrap-doc:Stores information about an isotopic cluster (i.e. potential peptide charge variants)
        IsotopeCluster(IsotopeCluster &) nogil except + # compiler

        ChargedIndexSet peaks
        libcpp_vector[ size_t ] scans

cdef extern from "<OpenMS/DATASTRUCTURES/IsotopeCluster.h>" namespace "OpenMS::IsotopeCluster":
    
    cdef cppclass ChargedIndexSet "OpenMS::IsotopeCluster::ChargedIndexSet":
        ChargedIndexSet() nogil except + # wrap-doc:Index set with associated charge estimate
        ChargedIndexSet(ChargedIndexSet) nogil except + #wrap-ignore
        Int charge

