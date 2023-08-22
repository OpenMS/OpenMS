from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>" namespace "OpenMS":
    
    cdef cppclass IntensityBalanceFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        IntensityBalanceFilter() except + nogil  # wrap-doc:It divides the m/z-range into ten regions and sums the intensity in these region
        IntensityBalanceFilter(IntensityBalanceFilter &) except + nogil 

        double apply(MSSpectrum & ) except + nogil 
        # POINTER # FilterFunctor * create() except + nogil 
        String getProductName() except + nogil 

