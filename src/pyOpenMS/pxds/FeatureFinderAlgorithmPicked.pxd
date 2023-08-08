from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmPicked(DefaultParamHandler):

        # wrap-inherits:
        #   DefaultParamHandler
        FeatureFinderAlgorithmPicked() except + nogil 
        # private
        FeatureFinderAlgorithmPicked(FeatureFinderAlgorithmPicked &) except + nogil  # wrap-ignore

        void setData(MSExperiment & input, FeatureMap & output, FeatureFinder & ff) except + nogil 
        void run() except + nogil 


        void setSeeds(FeatureMap& seeds) except + nogil 

        # static FeatureFinderAlgorithm* create()

#
# static methods are wrapped like this:
#

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>" namespace "OpenMS::FeatureFinderAlgorithmPicked":

    String getProductName()   except + nogil  # wrap-attach:FeatureFinderAlgorithmPicked
