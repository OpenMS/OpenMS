from ProgressLogger cimport *
from MSExperiment cimport *
from MSSpectrum cimport *

cdef extern from "<OpenMS/PROCESSING/SMOOTHING/ModifiedSincSmoother.h>" namespace "OpenMS":

    cdef cppclass ModifiedSincSmoother(ProgressLogger):
        # wrap-inherits:
        #   ProgressLogger

        ModifiedSincSmoother(bool isMS1, int degree, int m) except + nogil

        void filter(MSSpectrum & spectrum) except + nogil  # wrap-doc:Smoothes an MSSpectrum containing profile data
        void filterExperiment(MSExperiment & exp) except + nogil  # wrap-doc:Smoothes an MSExperiment containing profile data
