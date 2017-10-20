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
        #    DefaultParamHandler

        MarkerMower()            nogil except +
        MarkerMower(MarkerMower) nogil except + #wrap-ignore

        void filterSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakSpectrum(MSSpectrum & spec) nogil except +
        void filterPeakMap(MSExperiment & exp) nogil except +

        String getProductName() nogil except +

        void insertmarker(PeakMarker * peak_marker) nogil except +

