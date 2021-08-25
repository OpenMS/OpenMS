from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>" namespace "OpenMS":
    
    cdef cppclass ComplementFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor

        ComplementFilter() nogil except + # wrap-doc:Total intensity of peak pairs that could result from complementing fragments of charge state 1
        ComplementFilter(ComplementFilter &) nogil except +
        double apply(MSSpectrum & ) nogil except + # wrap-doc:Returns the total intensity of peak pairs which could result from complementing fragments

        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except + # wrap-doc:Returns the name for registration at the factory
