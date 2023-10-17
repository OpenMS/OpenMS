from Types cimport *
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>" namespace "OpenMS":
    
    cdef cppclass FilterFunctor(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        FilterFunctor() except + nogil   # wrap-doc:A FilterFunctor extracts some spectrum characteristics for quality assessment
        FilterFunctor(FilterFunctor &) except + nogil 

        # double apply(MSSpectrum & ) except + nogil 
        void registerChildren() except + nogil 
