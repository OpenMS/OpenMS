from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/ComplementFilter.h>" namespace "OpenMS":
    
    cdef cppclass ComplementFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor

        ComplementFilter() except + nogil  # wrap-doc:Total intensity of peak pairs that could result from complementing fragments of charge state 1
        ComplementFilter(ComplementFilter &) except + nogil 
        double apply(MSSpectrum & ) except + nogil  # wrap-doc:Returns the total intensity of peak pairs which could result from complementing fragments

        # POINTER # FilterFunctor * create() except + nogil 
        String getProductName() except + nogil  # wrap-doc:Returns the name for registration at the factory
