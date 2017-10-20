from Types cimport *

cdef extern from "<OpenMS/METADATA/Tagging.h>" namespace "OpenMS":
    
    cdef cppclass Tagging:

        Tagging() nogil except + 
        Tagging(Tagging) nogil except + 

        double getMassShift() nogil except +
        void setMassShift(double mass_shift) nogil except +

        IsotopeVariant getVariant() nogil except +
        void setVariant(IsotopeVariant variant) nogil except +

cdef extern from "<OpenMS/METADATA/Tagging.h>" namespace "OpenMS::Tagging":
    
    cdef enum IsotopeVariant "OpenMS::Tagging::IsotopeVariant":
        LIGHT
        HEAVY
        SIZE_OF_ISOTOPEVARIANT

