
from MSSpectrum cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>" namespace "OpenMS":

    cdef cppclass SpectrumAlignmentScore(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        SpectrumAlignmentScore() nogil except +
        SpectrumAlignmentScore(SpectrumAlignmentScore) nogil except +

        double operator()(MSSpectrum &, MSSpectrum &) nogil except + #wrap-ignore
        double operator()(MSSpectrum &) nogil except + #wrap-ignore
