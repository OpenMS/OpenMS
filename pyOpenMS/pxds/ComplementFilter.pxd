from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>" namespace "OpenMS":
    
    cdef cppclass ComplementFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        ComplementFilter() nogil except +
        ComplementFilter(ComplementFilter) nogil except +
        double apply(MSSpectrum[Peak1D] & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except +

