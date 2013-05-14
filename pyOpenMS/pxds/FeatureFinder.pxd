from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

from FeatureMap cimport *
from Feature cimport *
from Param cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>" namespace "OpenMS":

    cdef cppclass FeatureFinder(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger
        #

        FeatureFinder()      nogil except +
        void run(String algorithm_name, MSExperiment[Peak1D, ChromatogramPeak] & input_map,\
                 FeatureMap[Feature] & feats, Param & param, FeatureMap[Feature] & seeds) nogil except +

        Param getParameters(String algorithm_name)
