from libcpp cimport bool
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/ChromatogramPeak.h>" namespace "OpenMS::ChromatogramPeak":

    ctypedef double IntensityType
    ctypedef double CoordinateType
    ctypedef DPosition1 PositionType

cdef extern from "<OpenMS/KERNEL/ChromatogramPeak.h>" namespace "OpenMS":

    cdef cppclass ChromatogramPeak:

        ChromatogramPeak() nogil except +
        ChromatogramPeak(ChromatogramPeak &) nogil except +
        ChromatogramPeak(PositionType retention_time, IntensityType intensity) nogil except +
        bool operator==(ChromatogramPeak) nogil except +
        bool operator!=(ChromatogramPeak) nogil except +

        # We will not catch C++ exceptions for get/set methods for performance
        # reasons (no memory allocation is involved).

        IntensityType getIntensity() nogil 
        void setIntensity(IntensityType) nogil 

        DPosition1 getPosition() nogil 
        void setPosition(DPosition1) nogil 

        CoordinateType getRT() nogil 
        void setRT(CoordinateType) nogil 

        CoordinateType getPos() nogil 
        void setPos(CoordinateType) nogil 

        # alias for getRT 
        CoordinateType getMZ() nogil 
        void setMZ(CoordinateType) nogil 

