from TransitionTSVFile cimport *
from TargetedExperiment cimport *
from DefaultParamHandler cimport *
from GaussFilter cimport *
from SavitzkyGolayFilter cimport *
from MzMLFile cimport *
from MSExperiment cimport *
from FeatureMap cimport *
from PeakPickerHiRes cimport *

from String cimport *
from libcpp.vector cimport vector as libcpp_vector

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>" namespace "OpenMS":

    cdef cppclass TargetedSpectraExtractor(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        TargetedSpectraExtractor() nogil except +
        TargetedSpectraExtractor(TargetedSpectraExtractor) nogil except +

        void getDefaultParameters(Param) nogil except +

        void annotateSpectra(libcpp_vector[ MSSpectrum ], TargetedExperiment, libcpp_vector[ MSSpectrum ], FeatureMap) nogil except +

        void pickSpectrum(MSSpectrum, MSSpectrum) nogil except +

        void scoreSpectra(libcpp_vector[ MSSpectrum ], libcpp_vector[ MSSpectrum ], FeatureMap, libcpp_vector[ MSSpectrum ]) nogil except +

        void selectSpectra(libcpp_vector[ MSSpectrum ], FeatureMap, libcpp_vector[ MSSpectrum ], FeatureMap) nogil except +

        void selectSpectra(libcpp_vector[ MSSpectrum ], libcpp_vector[ MSSpectrum ]) nogil except +

        void extractSpectra(MSExperiment, TargetedExperiment, libcpp_vector[ MSSpectrum ], FeatureMap) nogil except +
