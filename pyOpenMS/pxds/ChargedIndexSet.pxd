from Types cimport *

cdef extern from "<OpenMS/DATASTRUCTURES/IsotopeCluster.h>" namespace "OpenMS::IsotopeCluster":
    
    cdef cppclass ChargedIndexSet "OpenMS::IsotopeCluster::ChargedIndexSet":
        ChargedIndexSet() nogil except +
        ChargedIndexSet(ChargedIndexSet) nogil except + #wrap-ignore
        Int charge

