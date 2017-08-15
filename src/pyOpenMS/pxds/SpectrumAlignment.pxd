from Types cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumAlignment(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        SpectrumAlignment() nogil except +
        SpectrumAlignment(SpectrumAlignment) nogil except +

        void getSpectrumAlignment(libcpp_vector[ libcpp_pair[ Size, Size ] ] & alignment, MSSpectrum & s1, MSSpectrum & s2) nogil except +  # wrap-ignore
