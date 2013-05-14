from libcpp cimport *
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/Peak2D.h>" namespace "OpenMS":

    cdef cppclass Peak2D:
        Peak2D()               nogil except +
        Peak2D(Peak2D &)               nogil except +
        Real getIntensity()     nogil except +
        DoubleReal getMZ()     nogil except +
        DoubleReal getRT()     nogil except +
        void setMZ(DoubleReal)  nogil except +
        void setRT(DoubleReal)  nogil except +
        void setIntensity(Real) nogil except +
        bool operator==(Peak2D) nogil except +
        bool operator!=(Peak2D) nogil except +
        # DPosition2  getPosition()
        # void setPosition(DPosition2 & position)

