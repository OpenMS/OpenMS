from libcpp cimport *
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/Peak1D.h>" namespace "OpenMS":

    cdef cppclass Peak1D:
        Peak1D() except + nogil 
        Peak1D(Peak1D &) except + nogil 

        # We will not catch C++ exceptions for get/set methods for performance
        # reasons (no memory allocation is involved).
        float getIntensity() nogil 
        double getMZ() nogil 
        void setMZ(double) nogil 
        void setIntensity(float) nogil 
        bool operator==(Peak1D) except + nogil 
        bool operator!=(Peak1D) except + nogil 
        double getPos() nogil 
        void setPos(double pos) nogil 
        # DPosition1 getPosition() except + nogil  # wrap-ignore
        # void setPosition(DPosition1 position) except + nogil  # wrap-ignore
    
