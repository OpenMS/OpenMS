from Types cimport *
from libcpp cimport bool
from MSSpectrum cimport *
from Peak1D cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/HyperScore.h>" namespace "OpenMS":
    
    cdef cppclass HyperScore "OpenMS::HyperScore":

        # compiler
        HyperScore() nogil except + # wrap-doc:An implementation of the X!Tandem HyperScore PSM scoring function
        HyperScore(HyperScore &) nogil except + # compiler

        double compute(double fragment_mass_tolerance, 
                       bool fragment_mass_tolerance_unit_ppm,
                       MSSpectrum & exp_spectrum, MSSpectrum & theo_spectrum) nogil except +
            # wrap-doc:
            #   Compute the (ln transformed) X!Tandem HyperScore 
            #   -----
            #   1. the dot product of peak intensities between matching peaks in experimental and theoretical spectrum is calculated
            #   2. the HyperScore is calculated from the dot product by multiplying by factorials of matching b- and y-ions
            #   -----
            #   :note: Peak intensities of the theoretical spectrum are typically 1 or TIC normalized, but can also be e.g. ion probabilities
            #   :param fragment_mass_tolerance: Mass tolerance applied left and right of the theoretical spectrum peak position
            #   :param fragment_mass_tolerance_unit_ppm: Unit of the mass tolerance is: Thomson if false, ppm if true
            #   :param exp_spectrum: Measured spectrum
            #   :param theo_spectrum: Theoretical spectrum Peaks need to contain an ion annotation as provided by TheoreticalSpectrumGenerator
