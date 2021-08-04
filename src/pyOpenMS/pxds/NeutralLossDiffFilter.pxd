from Types cimport *
from FilterFunctor cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass NeutralLossDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        # wrap-doc:
        #  NeutralLossDiffFilter returns the total intensity ob peak pairs whose m/z difference can be explained by a neutral loss

        NeutralLossDiffFilter() nogil except +
        NeutralLossDiffFilter(NeutralLossDiffFilter &) nogil except +
        double apply(MSSpectrum & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except +

