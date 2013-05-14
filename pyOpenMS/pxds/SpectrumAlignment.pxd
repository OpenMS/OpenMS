from Types cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumAlignment(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        SpectrumAlignment() nogil except +
        SpectrumAlignment(SpectrumAlignment) nogil except +
        # TODO nested STL
        # TEMPLATE # void getSpectrumAlignment(libcpp_vector[ libcpp_pair[ Size, Size ] ] & alignment, SpectrumType & s1, SpectrumType & s2) nogil except +

