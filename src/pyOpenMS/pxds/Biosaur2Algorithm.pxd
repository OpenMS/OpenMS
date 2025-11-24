from Types cimport *
from DefaultParamHandler cimport *
from FeatureMap cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>" namespace "OpenMS":

    cdef cppclass Biosaur2Algorithm(DefaultParamHandler):
        # wrap-doc:
        #  C++ implementation of the Biosaur2 feature detection workflow.

        Biosaur2Algorithm() except + nogil  # wrap-doc:Instantiate the Biosaur2 feature finding algorithm
        Biosaur2Algorithm(Biosaur2Algorithm&) except + nogil  # wrap-doc:Copy constructor

        void run(const MSExperiment& input,
                 FeatureMap& feature_map) except + nogil  # wrap-doc:Run the algorithm storing only the resulting features
