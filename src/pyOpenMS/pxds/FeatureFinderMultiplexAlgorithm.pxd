from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from ConsensusMap cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FEATUREFINDER/FeatureFinderMultiplexAlgorithm.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderMultiplexAlgorithm(DefaultParamHandler):

        # wrap-inherits:
        #   DefaultParamHandler
        FeatureFinderMultiplexAlgorithm() except + nogil 

        FeatureFinderMultiplexAlgorithm(FeatureFinderMultiplexAlgorithm &) except + nogil  # compiler

        void run(MSExperiment& exp, bool progress) except + nogil  # wrap-doc:Main method for feature detection

        FeatureMap getFeatureMap() except + nogil  # TODO

        ConsensusMap getConsensusMap() except + nogil  # TODO
