from Types cimport *
from OpenSwathDataStructures cimport *
from SpectrumAccessTransforming cimport *
from SpectrumAccessOpenMS cimport *
from SpectrumAccessOpenMSCached cimport *
from SpectrumAccessOpenMSInMemory cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumAccessQuadMZTransforming(SpectrumAccessTransforming) :
        # wrap-inherits:
        #  SpectrumAccessTransforming

        SpectrumAccessQuadMZTransforming() except + nogil  # wrap-pass-constructor
        SpectrumAccessQuadMZTransforming(SpectrumAccessQuadMZTransforming &) except + nogil  # compiler

        # SpectrumAccessQuadMZTransforming(shared_ptr[ ISpectrumAccess] sptr, double a, double b, double c, bool ppm) except + nogil 
        SpectrumAccessQuadMZTransforming(shared_ptr[ SpectrumAccessOpenMS ], double a, double b, double c, bool ppm) except + nogil 
        SpectrumAccessQuadMZTransforming(shared_ptr[ SpectrumAccessOpenMSCached ], double a, double b, double c, bool ppm) except + nogil 
        SpectrumAccessQuadMZTransforming(shared_ptr[ SpectrumAccessOpenMSInMemory ], double a, double b, double c, bool ppm) except + nogil 

