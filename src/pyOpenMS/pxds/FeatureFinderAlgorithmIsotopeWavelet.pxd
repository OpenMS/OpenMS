from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmIsotopeWavelet(DefaultParamHandler):

        # wrap-inherits:
        #    DefaultParamHandler
        FeatureFinderAlgorithmIsotopeWavelet() nogil except +

        void setData(MSExperiment & input, FeatureMap& output, FeatureFinder & ff) nogil except +
        void run() nogil except +
        # MSSpectrum * createHRData(const UInt i)
        # static FeatureFinderAlgorithm* create()

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>" namespace "OpenMS::FeatureFinderAlgorithmIsotopeWavelet":

    String getProductName()   nogil except + # wrap-attach:FeatureFinderAlgorithmIsotopeWavelet
