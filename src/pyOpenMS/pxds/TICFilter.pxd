from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>" namespace "OpenMS":
    
    cdef cppclass TICFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        TICFilter() nogil except +
        TICFilter(TICFilter) nogil except +
        double apply(MSSpectrum[Peak1D] & ) nogil except +
        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except +

