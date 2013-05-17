from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        IsotopeDiffFilter() nogil except +
        IsotopeDiffFilter(IsotopeDiffFilter) nogil except +
        double apply(MSSpectrum[Peak1D] & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except +

