from libcpp cimport bool
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/ChromatogramPeak.h>" namespace "OpenMS::ChromatogramPeak":

    ctypedef double IntensityType
    ctypedef double CoordinateType
    ctypedef DPosition1 PositionType

cdef extern from "<OpenMS/KERNEL/ChromatogramPeak.h>" namespace "OpenMS":

    cdef cppclass ChromatogramPeak:

        ChromatogramPeak() nogil except + # wrap-doc:A 1-dimensional raw data point or peak for chromatograms
        ChromatogramPeak(ChromatogramPeak &) nogil except +
        ChromatogramPeak(PositionType retention_time, IntensityType intensity) nogil except +
        bool operator==(ChromatogramPeak) nogil except +
        bool operator!=(ChromatogramPeak) nogil except +

        # We will not catch C++ exceptions for get/set methods for performance
        # reasons (no memory allocation is involved).

        IntensityType getIntensity() nogil except + # wrap-doc:Non-mutable access to the data point intensity (height)
        void setIntensity(IntensityType) nogil except + # wrap-doc:Mutable access to the data point intensity (height)

        DPosition1 getPosition() nogil except + # TODO
        void setPosition(DPosition1) nogil except + # TODO

        CoordinateType getRT() nogil except + # wrap-doc:Non-mutable access to RT
        void setRT(CoordinateType) nogil except + # wrap-doc:Mutable access to RT

        CoordinateType getPos() nogil except + # wrap-doc:Alias for getRT()
        void setPos(CoordinateType) nogil except + # wrap-doc:Alias for setRT()

        # alias for getRT 
        CoordinateType getMZ() nogil except + # wrap-doc:Alias for getRT()
        void setMZ(CoordinateType) nogil except + # wrap-doc:Alias for setRT()

