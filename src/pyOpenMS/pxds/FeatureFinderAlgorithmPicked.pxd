from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmPicked(DefaultParamHandler):
        # wrap-doc:
        #  The purpose of this struct is to provide definitions of classes and
        #  typedefs which are used throughout all FeatureFinder classes. */
        #  struct OPENMS_DLLAPI FeatureFinderDefs { /// Index to peak consisting
        #  of two UInts (scan index / peak index) typedef
        #  IsotopeCluster::IndexPair IndexPair;

        # wrap-inherits:
        #   DefaultParamHandler
        FeatureFinderAlgorithmPicked() except + nogil 
        # private
        FeatureFinderAlgorithmPicked(FeatureFinderAlgorithmPicked &) except + nogil  # wrap-ignore

        void run(MSExperiment & input_map, FeatureMap & output, Param & param, FeatureMap & seeds) except + nogil



