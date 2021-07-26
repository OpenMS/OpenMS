from libcpp cimport bool
from Types cimport *
from Peak2D cimport *
from MetaInfoInterface cimport *
from UniqueIdInterface cimport *

cdef extern from "<OpenMS/KERNEL/RichPeak2D.h>" namespace "OpenMS":


    cdef cppclass RichPeak2D(Peak2D, UniqueIdInterface, MetaInfoInterface):
        # wrap-inherits:
        #    Peak2D
        #    UniqueIdInterface
        #    MetaInfoInterface

        RichPeak2D() nogil except +
        RichPeak2D(RichPeak2D &) nogil except +
        #RichPeak2D(DPosition2 &, float) nogil except +

        bool operator==(RichPeak2D) nogil except +
        bool operator!=(RichPeak2D) nogil except +

