from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmSH[PeakT, FeatureT](DefaultParamHandler):

        # wrap-inherits:
        #    DefaultParamHandler
        #
        # wrap-instances:
        #    FeatureFinderAlgorithmSH := FeatureFinderAlgorithmSH[Peak1D, Feature]

        FeatureFinderAlgorithmSH()      nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>" namespace "OpenMS::FeatureFinderAlgorithmSH<OpenMS::Peak1D,OpenMS::Feature>":

    String getProductName()   nogil except + # wrap-attach:FeatureFinderAlgorithmSH
