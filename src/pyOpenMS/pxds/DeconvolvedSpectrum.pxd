from MSSpectrum cimport *
from FLASHHelperClasses cimport *



cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>" namespace "OpenMS":

    cdef cppclass DeconvolvedSpectrum:
        # wrap-inherits:

        # default constructor
        DeconvolvedSpectrum() except + nogil
        # copy constructor
        DeconvolvedSpectrum(DeconvolvedSpectrum &) except + nogil

