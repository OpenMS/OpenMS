from libcpp cimport *
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/Peak1D.h>" namespace "OpenMS":

    cdef cppclass Peak1D:
        Peak1D() nogil except +
        Peak1D(Peak1D &) nogil except +

        # We will not catch C++ exceptions for get/set methods for performance
        # reasons (no memory allocation is involved).
        float getIntensity() nogil 
        double getMZ() nogil 
        void setMZ(double) nogil 
        void setIntensity(float) nogil 
        bool operator==(Peak1D) nogil except +
        bool operator!=(Peak1D) nogil except +
        double getPos() nogil 
        void setPos(double pos) nogil 
        # DPosition1 getPosition() nogil except + # wrap-ignore
        # void setPosition(DPosition1 position) nogil except + # wrap-ignore
    
