from libcpp cimport bool
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/ChromatogramPeak.h>" namespace "OpenMS::ChromatogramPeak":

    ctypedef DoubleReal IntensityType
    ctypedef DoubleReal CoordinateType

cdef extern from "<OpenMS/KERNEL/ChromatogramPeak.h>" namespace "OpenMS":

    cdef cppclass ChromatogramPeak:

        ChromatogramPeak()     nogil except +
        ChromatogramPeak(ChromatogramPeak) nogil except +
        bool operator==(ChromatogramPeak) nogil except +
        bool operator!=(ChromatogramPeak) nogil except +

        IntensityType getIntensity() nogil except +
        void setIntensity(IntensityType) nogil except +

        #DPosition1 PositionType getPosition() nogil except +
        void setPosition(DPosition1) nogil except +

        CoordinateType getRT() nogil except +
        void setRT(CoordinateType) nogil except +

        CoordinateType getPos() nogil except +
        void setPos(CoordinateType) nogil except +

        # alias for getRT 
        CoordinateType getMZ() nogil except +
        void setMZ(CoordinateType) nogil except +
