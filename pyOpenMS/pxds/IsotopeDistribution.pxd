from libcpp cimport bool
from Types cimport *
from String cimport *

cdef extern from "<OpenMS/CHEMISTRY/IsotopeDistribution.h>" namespace "OpenMS":

    cdef cppclass IsotopeDistribution:

        IsotopeDistribution() nogil except +
        IsotopeDistribution(IsotopeDistribution) nogil except + # wrap-ignore

