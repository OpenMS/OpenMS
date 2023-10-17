from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>" namespace "OpenMS":
    
    cdef cppclass TICFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        TICFilter() except + nogil  # TODO
        TICFilter(TICFilter &) except + nogil 
        double apply(MSSpectrum & ) except + nogil 
        # POINTER # FilterFunctor * create() except + nogil 
        String getProductName() except + nogil 
