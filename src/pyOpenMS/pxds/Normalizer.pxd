# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/PROCESSING/SCALING/Normalizer.h>" namespace "OpenMS":

    cdef cppclass Normalizer(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler
        # wrap-doc:
        #  Normalizes the peak intensities spectrum-wise

        Normalizer() except + nogil 

        Normalizer(Normalizer) except + nogil 

        void filterSpectrum(MSSpectrum & spec) except + nogil  # wrap-doc:Normalizes the spectrum

        void filterPeakSpectrum(MSSpectrum & spec) except + nogil  # wrap-doc:Normalizes the peak spectrum

        void filterPeakMap(MSExperiment & exp) except + nogil  # wrap-doc:Normalizes the peak map
