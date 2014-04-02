from libcpp cimport *
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/Peak1D.h>" namespace "OpenMS":

    cdef cppclass Peak1D:
        Peak1D()               nogil except +
        Peak1D(Peak1D &)               nogil except +
        float getIntensity()     nogil except +
        double getMZ()     nogil except +
        void setMZ(double)  nogil except +
        void setIntensity(float) nogil except +
        bool operator==(Peak1D) nogil except +
        bool operator!=(Peak1D) nogil except +
        double getPos() nogil except +
        void setPos(double pos) nogil except +
        # DPosition1 getPosition() nogil except + # wrap-ignore
        # void setPosition(DPosition1 position) nogil except + # wrap-ignore
    
