from Types cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>" namespace "OpenMS":
    
    cdef cppclass FilterFunctor(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        FilterFunctor() nogil except +
        FilterFunctor(FilterFunctor) nogil except +
        # double apply(MSSpectrum & ) nogil except +
        void registerChildren() nogil except +

