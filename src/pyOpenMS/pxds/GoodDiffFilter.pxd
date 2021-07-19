from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass GoodDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        GoodDiffFilter() nogil except + # wrap-doc:It counts the number ob peak pairs whose m/z difference can be explained by a amino acid loss
        GoodDiffFilter(GoodDiffFilter) nogil except + # wrap-ignore
        double apply(MSSpectrum & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except + # TODO

