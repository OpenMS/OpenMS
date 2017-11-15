from TransitionTSVReader cimport *
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

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/SpectrumExtractor.h>" namespace "OpenMS":

    cdef cppclass SpectrumExtractor(DefaultParamHandler):
        # wrap-inherits:
        #  DefaultParamHandler

        SpectrumExtractor() nogil except +
        SpectrumExtractor(SpectrumExtractor) nogil except +

        void setRTWindow(double) nogil except +
        double getRTWindow() nogil except +

        void setMinScore(double) nogil except +
        double getMinScore() nogil except +

        void setMinForwardMatch(double) nogil except +
        double getMinForwardMatch() nogil except +

        void setMinReverseMatch(double) nogil except +
        double getMinReverseMatch() nogil except +

        void setMZTolerance(double) nogil except +
        double getMZTolerance() nogil except +

        void setMZToleranceUnits(String) nogil except +
        String getMZToleranceUnits() nogil except +

        void setSGolayFrameLength(UInt) nogil except +
        UInt getSGolayFrameLength() nogil except +

        void setSGolayPolynomialOrder(UInt) nogil except +
        UInt getSGolayPolynomialOrder() nogil except +

        void setGaussWidth(double) nogil except +
        double getGaussWidth() nogil except +

        void setUseGauss(bool) nogil except +
        bool getUseGauss() nogil except +

        void setSignalToNoise(double) nogil except +
        double getSignalToNoise() nogil except +

        void setPeakHeightMin(double) nogil except +
        double getPeakHeightMin() nogil except +

        void setPeakHeightMax(double) nogil except +
        double getPeakHeightMax() nogil except +

        void setFWHMThreshold(double) nogil except +
        double getFWHMThreshold() nogil except +

        void setTICWeight(double) nogil except +
        double getTICWeight() nogil except +

        void setFWHMWeight(double) nogil except +
        double getFWHMWeight() nogil except +

        void setSNRWeight(double) nogil except +
        double getSNRWeight() nogil except +

        void getDefaultParameters(Param) nogil except +

        void annotateSpectra(libcpp_vector[ MSSpectrum ], TargetedExperiment, libcpp_vector[ MSSpectrum ], FeatureMap) nogil except +

        void pickSpectrum(MSSpectrum, MSSpectrum) nogil except +

        void scoreSpectra(libcpp_vector[ MSSpectrum ], libcpp_vector[ MSSpectrum ], FeatureMap, libcpp_vector[ MSSpectrum ]) nogil except +

        void selectSpectra(libcpp_vector[ MSSpectrum ], FeatureMap, libcpp_vector[ MSSpectrum ], FeatureMap) nogil except +

        void selectSpectra(libcpp_vector[ MSSpectrum ], libcpp_vector[ MSSpectrum ]) nogil except +

        void extractSpectra(MSExperiment, TargetedExperiment, libcpp_vector[ MSSpectrum ], FeatureMap) nogil except +
