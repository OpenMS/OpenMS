from Types cimport *
from FilterFunctor cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/IsotopeDiffFilter.h>" namespace "OpenMS":
    
    cdef cppclass IsotopeDiffFilter(FilterFunctor) :
        # wrap-inherits:
        #  FilterFunctor
        IsotopeDiffFilter() nogil except + # wrap-doc:IsotopeDiffFilter returns total intensity of peak pairs that could result from isotope peaks
        IsotopeDiffFilter(IsotopeDiffFilter &) nogil except +
        double apply(MSSpectrum & ) nogil except + # TODO

        # POINTER # FilterFunctor * create() nogil except +
        String getProductName() nogil except + # TODO

