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

        IntensityType getIntensity() nogil except +
        void setIntensity(IntensityType) nogil except +

        DPosition1 getPosition() nogil except +
        void setPosition(DPosition1) nogil except +

        CoordinateType getRT() nogil except +
        void setRT(CoordinateType) nogil except +

        CoordinateType getPos() nogil except +
        void setPos(CoordinateType) nogil except +

        # alias for getRT 
        CoordinateType getMZ() nogil except +
        void setMZ(CoordinateType) nogil except +
