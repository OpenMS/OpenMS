# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *
from PeakMarker cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/MarkerMower.h>" namespace "OpenMS":

    cdef cppclass MarkerMower(DefaultParamHandler):
        # wrap-inherits:
        #   DefaultParamHandler

        MarkerMower() except + nogil 
        MarkerMower(MarkerMower &) except + nogil 

        void filterSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakSpectrum(MSSpectrum & spec) except + nogil 
        void filterPeakMap(MSExperiment & exp) except + nogil 

        String getProductName() except + nogil  # wrap-doc:Returns the product name

        void insertmarker(PeakMarker * peak_marker) except + nogil  # wrap-doc:Insert new Marker (violates the DefaultParamHandler interface)
