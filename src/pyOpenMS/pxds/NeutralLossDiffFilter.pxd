from Types cimport *
from FilterFunctor cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass NeutralLossDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        # wrap-doc:
        #  NeutralLossDiffFilter returns the total intensity ob peak pairs whose m/z difference can be explained by a neutral loss

        NeutralLossDiffFilter() except + nogil 
        NeutralLossDiffFilter(NeutralLossDiffFilter &) except + nogil 
        double apply(MSSpectrum & ) except + nogil 
        # POINTER # FilterFunctor * create() except + nogil 
        String getProductName() except + nogil 

