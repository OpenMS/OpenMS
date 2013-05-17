from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>" namespace "OpenMS":
    
    cdef cppclass IntensityBalanceFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        IntensityBalanceFilter() nogil except +
        IntensityBalanceFilter(IntensityBalanceFilter) nogil except +
        double apply(MSSpectrum[Peak1D] & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except +

