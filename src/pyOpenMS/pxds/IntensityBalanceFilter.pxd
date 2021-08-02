from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>" namespace "OpenMS":
    
    cdef cppclass IntensityBalanceFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        IntensityBalanceFilter() nogil except + # wrap-doc:It divides the m/z-range into ten regions and sums the intensity in these region
        IntensityBalanceFilter(IntensityBalanceFilter &) nogil except +

        double apply(MSSpectrum & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except +

