from libcpp cimport *
from Types cimport *
from DPosition cimport *

cdef extern from "<OpenMS/KERNEL/MobilityPeak1D.h>" namespace "OpenMS":

    cdef cppclass MobilityPeak1D:
        # wrap-hash:
        #  std
        MobilityPeak1D() except + nogil
        MobilityPeak1D(MobilityPeak1D &) except + nogil 

        # We will not catch C++ exceptions for get/set methods for performance
        # reasons (no memory allocation is involved).
        float getIntensity() nogil 
        double getMobility() nogil 
        void setMobility(double) nogil 
        void setIntensity(float) nogil 
        bool operator==(MobilityPeak1D) except + nogil 
        bool operator!=(MobilityPeak1D) except + nogil 
        double getPos() nogil 
        void setPos(double pos) nogil 
        # DPosition1 getPosition() except + nogil  # wrap-ignore
        # void setPosition(DPosition1 position) except + nogil  # wrap-ignore
    
