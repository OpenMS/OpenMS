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

        FeatureFinder() nogil except +

        FeatureFinder(FeatureFinder &) nogil except + # compiler

        void run(String algorithm_name,
                 MSExperiment & input_map,
                 FeatureMap & feats,
                 Param & param,
                 FeatureMap & seeds
                 ) nogil except +
            # wrap-doc:
            #   Executes the FeatureFinder using the given algorithm
            #   -----
            #   There are several constraints for the `input_map`.  They are tested before the algorithm starts.  It must only contain MS 1 level scans and you have to call updateRanges() before passing it to this method
            #   The input map is sorted by RT & m/z if that's not the case
            #   Furthermore we throw an Exception if the data contains negative m/z values,
            #   as this will disturb most algorithms
            #   -----
            #   :param algorithm_name: Name of the feature finding algorithm to use
            #   :param input_map: Input peak map
            #   :param features: Output feature map
            #   :param param: Algorithm parameters
            #   :param seeds: List of seeds to use

        Param getParameters(String algorithm_name) nogil except + # wrap-doc:Returns the default parameters for the algorithm with name `algorithm_name`
