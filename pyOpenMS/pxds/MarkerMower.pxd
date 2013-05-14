# part of the SpectraFilters
from MSSpectrum cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FILTERING/TRANSFORMERS/MarkerMower.h>" namespace "OpenMS":

    cdef cppclass MarkerMower(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        MarkerMower()            nogil except +
        MarkerMower(MarkerMower) nogil except + #wrap-ignore

        void filterSpectrum(MSSpectrum[Peak1D] & spec) nogil except +
        void filterPeakSpectrum(MSSpectrum[Peak1D] & spec) nogil except +
        void filterPeakMap(MSExperiment[Peak1D, ChromatogramPeak] & exp) nogil except +

        String getProductName() nogil except +

        # Pointer handling, still TODO 
        # void insertmarker(PeakMarker * peak_marker);

