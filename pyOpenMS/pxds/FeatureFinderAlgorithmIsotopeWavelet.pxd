from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmIsotopeWavelet[PeakT, FeatureT](DefaultParamHandler):

        # wrap-inherits:
        #    DefaultParamHandler
        #
        # wrap-instances:
        #    FeatureFinderAlgorithmIsotopeWavelet := FeatureFinderAlgorithmIsotopeWavelet[Peak1D, Feature]

        FeatureFinderAlgorithmIsotopeWavelet()      nogil except +

        void setData(MSExperiment[Peak1D, ChromatogramPeak] & input, FeatureMap[Feature] & output, FeatureFinder & ff)
        void run()

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>" namespace "OpenMS::FeatureFinderAlgorithmIsotopeWavelet<OpenMS::Peak1D,OpenMS::Feature>":

    String getProductName()   nogil except + # wrap-attach:FeatureFinderAlgorithmIsotopeWavelet
