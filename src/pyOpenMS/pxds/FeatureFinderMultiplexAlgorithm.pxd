from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from ConsensusMap cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderMultiplexAlgorithm(DefaultParamHandler):

        # wrap-inherits:
        #    DefaultParamHandler
        FeatureFinderMultiplexAlgorithm()      nogil except +

        void run(MSExperiment& exp, bool progress) nogil except +

        FeatureMap getFeatureMap() nogil except +

        ConsensusMap getConsensusMap() nogil except +
