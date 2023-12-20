from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        IsotopeDiffFilter() except + nogil  # wrap-doc:IsotopeDiffFilter returns total intensity of peak pairs that could result from isotope peaks
        IsotopeDiffFilter(IsotopeDiffFilter &) except + nogil 
        double apply(MSSpectrum & ) except + nogil  # TODO

        # POINTER # FilterFunctor * create() except + nogil 
        String getProductName() except + nogil  # TODO

