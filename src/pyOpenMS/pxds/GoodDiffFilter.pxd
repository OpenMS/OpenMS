from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/GoodDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass GoodDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        GoodDiffFilter() nogil except +
        GoodDiffFilter(GoodDiffFilter) nogil except +
        double apply(MSSpectrum[Peak1D] & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except +

