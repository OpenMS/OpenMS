from Types cimport *
from DPosition cimport *
from libcpp cimport bool

cdef extern from "<OpenMS/KERNEL/RangeManager.h>" namespace "OpenMS":
    
    cdef cppclass RangeManager1 "OpenMS::RangeManager<1>":
        # wrap-ignore
        RangeManager1() nogil except +
        RangeManager1(RangeManager1) nogil except +
        DPosition1 getMin() nogil except +
        DPosition1 getMax() nogil except +
        double getMinInt() nogil except +
        double getMaxInt() nogil except +
        void clearRanges() nogil except +

    cdef cppclass RangeManager2 "OpenMS::RangeManager<2>":
        # wrap-ignore
        RangeManager2() nogil except +
        RangeManager2(RangeManager2) nogil except +
        DPosition2 getMin() nogil except +
        DPosition2 getMax() nogil except +
        double getMinInt() nogil except +
        double getMaxInt() nogil except +
        void clearRanges() nogil except +
