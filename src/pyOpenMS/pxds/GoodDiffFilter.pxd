from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass GoodDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        # wrap-doc:
        #  It counts the number ob peak pairs whose m/z difference can be explained by a amino acid loss

        GoodDiffFilter() except + nogil 
        GoodDiffFilter(GoodDiffFilter &) except + nogil 

        double apply(MSSpectrum & ) except + nogil 
        # POINTER # FilterFunctor * create() except + nogil 
        String getProductName() except + nogil  # wrap-doc:Returns the final product name 
