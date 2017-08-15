from Types cimport *
from libcpp cimport bool
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/HyperScore.h>" namespace "OpenMS":
    
    cdef cppclass HyperScore "OpenMS::HyperScore":
        HyperScore() nogil except + 
        HyperScore(HyperScore) nogil except + #wrap-ignore

        double compute(double fragment_mass_tolerance, 
                       bool fragment_mass_tolerance_unit_ppm,
                       MSSpectrum & exp_spectrum, MSSpectrum & theo_spectrum) nogil except +

