from libcpp cimport *
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/Peak2D.h>" namespace "OpenMS":

    cdef cppclass Peak2D:
        Peak2D()               nogil except +
        Peak2D(Peak2D &)               nogil except +
        float getIntensity()     nogil except +
        double getMZ()     nogil except +
        double getRT()     nogil except +
        void setMZ(double)  nogil except +
        void setRT(double)  nogil except +
        void setIntensity(float) nogil except +
        bool operator==(Peak2D) nogil except +
        bool operator!=(Peak2D) nogil except +
        # DPosition2  getPosition()
        # void setPosition(DPosition2 & position)

cdef extern from "<OpenMS/KERNEL/Peak2D.h>" namespace "OpenMS::Peak2D":
    
    cdef enum DimensionDescription "OpenMS::Peak2D::DimensionDescription":
        RT
        MZ
        DIMENSION
