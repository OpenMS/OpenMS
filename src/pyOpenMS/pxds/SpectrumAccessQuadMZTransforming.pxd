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

        SpectrumAccessQuadMZTransforming() nogil except + # wrap-pass-constructor
        SpectrumAccessQuadMZTransforming(SpectrumAccessQuadMZTransforming &) nogil except + # compiler

        # SpectrumAccessQuadMZTransforming(shared_ptr[ ISpectrumAccess] sptr, double a, double b, double c, bool ppm) nogil except +
        SpectrumAccessQuadMZTransforming(shared_ptr[ SpectrumAccessOpenMS ], double a, double b, double c, bool ppm) nogil except +
        SpectrumAccessQuadMZTransforming(shared_ptr[ SpectrumAccessOpenMSCached ], double a, double b, double c, bool ppm) nogil except +
        SpectrumAccessQuadMZTransforming(shared_ptr[ SpectrumAccessOpenMSInMemory ], double a, double b, double c, bool ppm) nogil except +

