from libcpp.vector cimport vector as libcpp_vector
from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/IsotopeCluster.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeCluster "OpenMS::IsotopeCluster":

        IsotopeCluster() except + nogil  # wrap-doc:Stores information about an isotopic cluster (i.e. potential peptide charge variants)
        IsotopeCluster(IsotopeCluster &) except + nogil  # compiler

        ChargedIndexSet peaks
        libcpp_vector[ size_t ] scans

cdef extern from "<OpenMS/DATASTRUCTURES/IsotopeCluster.h>" namespace "OpenMS::IsotopeCluster":
    
    cdef cppclass ChargedIndexSet "OpenMS::IsotopeCluster::ChargedIndexSet":
        ChargedIndexSet() except + nogil  # wrap-doc:Index set with associated charge estimate
        ChargedIndexSet(ChargedIndexSet) except + nogil  #wrap-ignore
        Int charge

