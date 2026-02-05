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
        # wrap-doc:
        #  FeatureFinderMultiplexAlgorithm is a tool for the fully automated
        #  analysis of quantitative proteomics data. It detects pairs of isotopic
        #  envelopes with fixed m/z separation. It requires no prior sequence
        #  identification of the peptides and works on both profile or centroided
        #  spectra. In what follows we outline the algorithm

        # wrap-inherits:
        #   DefaultParamHandler
        FeatureFinderMultiplexAlgorithm() except + nogil 

        FeatureFinderMultiplexAlgorithm(FeatureFinderMultiplexAlgorithm &) except + nogil  # compiler

        void run(MSExperiment& exp, bool progress) except + nogil  # wrap-doc:Main method for feature detection

        FeatureMap getFeatureMap() except + nogil  # TODO

        ConsensusMap getConsensusMap() except + nogil  # TODO
