from libcpp cimport bool
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/ChromatogramPeak.h>" namespace "OpenMS::ChromatogramPeak":

    ctypedef double IntensityType
    ctypedef double CoordinateType
    ctypedef DPosition1 PositionType

cdef extern from "<OpenMS/KERNEL/ChromatogramPeak.h>" namespace "OpenMS":

    cdef cppclass ChromatogramPeak:

        ChromatogramPeak() except + nogil  # wrap-doc:A 1-dimensional raw data point or peak for chromatograms
        ChromatogramPeak(ChromatogramPeak &) except + nogil 
        ChromatogramPeak(PositionType retention_time, IntensityType intensity) except + nogil 
        bool operator==(ChromatogramPeak) except + nogil 
        bool operator!=(ChromatogramPeak) except + nogil 

        # We will not catch C++ exceptions for get/set methods for performance
        # reasons (no memory allocation is involved).

        IntensityType getIntensity() except + nogil  # wrap-doc:Returns the intensity
        void setIntensity(IntensityType) except + nogil  # wrap-doc:Sets the intensity

        DPosition1 getPosition() except + nogil  # TODO
        void setPosition(DPosition1) except + nogil  # TODO

        CoordinateType getRT() except + nogil  # wrap-doc:Returns the retention time
        void setRT(CoordinateType) except + nogil  # wrap-doc:Sets retention time

        CoordinateType getPos() except + nogil  # wrap-doc:Alias for getRT()
        void setPos(CoordinateType) except + nogil  # wrap-doc:Alias for setRT()

        # alias for getRT!!! 
        # Current not exposed since it does the same as getRT()a
        # and is probably used somewhere internally for filtering.
        # CoordinateType getMZ() except + nogil  # wrap-doc:Alias for getRT()
        # void setMZ(CoordinateType) except + nogil  # wrap-doc:Alias for setRT()
