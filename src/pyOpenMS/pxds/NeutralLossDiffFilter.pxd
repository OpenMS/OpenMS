from Types cimport *
from FilterFunctor cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/NeutralLossDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass NeutralLossDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        NeutralLossDiffFilter() nogil except +
        NeutralLossDiffFilter(NeutralLossDiffFilter) nogil except +
        double apply(MSSpectrum[Peak1D] & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except +

