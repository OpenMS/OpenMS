from libcpp cimport bool
from Types cimport *
from Peak2D cimport *
from MetaInfoInterface cimport *
from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/RichPeak2D.h>" namespace "OpenMS":


    cdef cppclass RichPeak2D(Peak2D, UniqueIdInterface, MetaInfoInterface):
        # wrap-inherits:
        #   Peak2D
        #   UniqueIdInterface
        #   MetaInfoInterface

        RichPeak2D() except + nogil  # wrap-doc:A 2-dimensional raw data point or peak with meta information
        RichPeak2D(RichPeak2D &) except + nogil 
        #RichPeak2D(DPosition2 &, float) except + nogil 

        bool operator==(RichPeak2D) except + nogil 
        bool operator!=(RichPeak2D) except + nogil 

