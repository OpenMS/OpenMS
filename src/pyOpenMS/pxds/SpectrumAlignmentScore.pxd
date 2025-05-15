
from MSSpectrum cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/COMPARISON/SpectrumAlignmentScore.h>" namespace "OpenMS":

    cdef cppclass SpectrumAlignmentScore(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler
        SpectrumAlignmentScore() except + nogil 
        # wrap-doc:
                #  Similarity score via spectra alignment
                #  
                #  This class implements a simple scoring based on the alignment of spectra. This alignment
                #  is implemented in the SpectrumAlignment class and performs a dynamic programming alignment
                #  of the peaks, minimizing the distances between the aligned peaks and maximizing the number
                #  of peak pairs
                #  
                #  The scoring is done via the simple formula score = sum / (sqrt(sum1 * sum2)). sum is the
                #  product of the intensities of the aligned peaks, with the given exponent (default is 2)
                #  sum1 and sum2 are the sum of the intensities squared for each peak of both spectra respectively

        SpectrumAlignmentScore(SpectrumAlignmentScore &) except + nogil 

        double operator()(MSSpectrum &, MSSpectrum &) except + nogil  #wrap-ignore
        double operator()(MSSpectrum &) except + nogil  #wrap-ignore
