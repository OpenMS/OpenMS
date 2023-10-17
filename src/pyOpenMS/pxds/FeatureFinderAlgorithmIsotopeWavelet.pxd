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
        #   DefaultParamHandler
        FeatureFinderAlgorithmIsotopeWavelet() except + nogil 

        void setData(MSExperiment & input, FeatureMap& output, FeatureFinder & ff) except + nogil 
        void run() except + nogil 
        # MSSpectrum * createHRData(const UInt i)
        # static FeatureFinderAlgorithm* create()

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>" namespace "OpenMS::FeatureFinderAlgorithmIsotopeWavelet":

    String getProductName()   except + nogil  # wrap-attach:FeatureFinderAlgorithmIsotopeWavelet
