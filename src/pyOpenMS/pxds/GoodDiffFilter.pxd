from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass GoodDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        # wrap-doc:
        #   It counts the number ob peak pairs whose m/z difference can be explained by a amino acid loss

        GoodDiffFilter() nogil except +
        GoodDiffFilter(GoodDiffFilter &) nogil except +

        double apply(MSSpectrum & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except + # wrap-doc:Returns the final product name 
