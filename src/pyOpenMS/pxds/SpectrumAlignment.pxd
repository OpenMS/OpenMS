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
        # wrap-doc:
        #   Aligns the peaks of two sorted spectra
        #   Method 1: Using a banded (width via 'tolerance' parameter) alignment if absolute tolerances are given
        #       Scoring function is the m/z distance between peaks. Intensity does not play a role!
        #   Method 2: If relative tolerance (ppm) is specified a simple matching of peaks is performed:
        #   Peaks from s1 (usually the theoretical spectrum) are assigned to the closest peak in s2 if it lies in the tolerance window
        #   -----
        #   note: A peak in s2 can be matched to none, one or multiple peaks in s1. Peaks in s1 may be matched to none or one peak in s2
        #   note: Intensity is ignored 

        SpectrumAlignment() nogil except +
        SpectrumAlignment(SpectrumAlignment &) nogil except +

        void getSpectrumAlignment(libcpp_vector[ libcpp_pair[ Size, Size ] ] & alignment, MSSpectrum & s1, MSSpectrum & s2) nogil except +  # wrap-ignore
