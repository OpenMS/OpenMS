from Types cimport *

cdef extern from "<OpenMS/METADATA/Tagging.h>" namespace "OpenMS":
    
    cdef cppclass Tagging:
        # wrap-doc:
                #  Meta information about tagging of a sample e.g. ICAT labeling
                #  
                #  Holds information about the mass difference between light and heavy tag
                #  All other relevant information is provided by Modification

        Tagging() except + nogil  
        Tagging(Tagging &) except + nogil 

        double getMassShift() except + nogil  # wrap-doc:Returns the mass difference between light and heavy variant (default is 0.0)
        void setMassShift(double mass_shift) except + nogil  # wrap-doc:Sets the mass difference between light and heavy variant

        IsotopeVariant getVariant() except + nogil  # wrap-doc:Returns the isotope variant of the tag (default is LIGHT)
        void setVariant(IsotopeVariant variant) except + nogil  # wrap-doc:Sets the isotope variant of the tag

cdef extern from "<OpenMS/METADATA/Tagging.h>" namespace "OpenMS::Tagging":
    
    cdef enum IsotopeVariant "OpenMS::Tagging::IsotopeVariant":
        LIGHT
        HEAVY
        SIZE_OF_ISOTOPEVARIANT

