from Types cimport *
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from RichPeak1D cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>" namespace "OpenMS":
    
    cdef cppclass SpectrumAlignment(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        SpectrumAlignment() nogil except +
        SpectrumAlignment(SpectrumAlignment) nogil except +

        # those two are wrapped manually, as nested stl type of first arg is not supported
        # by autowrap (yet). see addons/SpectrumAlignment.pyx
        void getSpectrumAlignment(libcpp_vector[ libcpp_pair[ Size, Size ] ] & alignment, MSSpectrum[Peak1D] & s1, MSSpectrum[Peak1D] & s2) nogil except +  # wrap-ignore
        void getSpectrumAlignment(libcpp_vector[ libcpp_pair[ Size, Size ] ] & alignment, MSSpectrum[RichPeak1D] & s1, MSSpectrum[RichPeak1D] & s2) nogil except +  # wrap-ignore

