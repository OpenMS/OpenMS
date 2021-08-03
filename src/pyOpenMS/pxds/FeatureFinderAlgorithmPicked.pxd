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
        #    DefaultParamHandler
        FeatureFinderAlgorithmPicked() nogil except +
        # private
        FeatureFinderAlgorithmPicked(FeatureFinderAlgorithmPicked &) nogil except + # wrap-ignore

        void setData(MSExperiment & input, FeatureMap & output, FeatureFinder & ff) nogil except +
        void run() nogil except +


        void setSeeds(FeatureMap& seeds) nogil except +

        # static FeatureFinderAlgorithm* create()

#
# static methods are wrapped like this:
#

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>" namespace "OpenMS::FeatureFinderAlgorithmPicked":

    String getProductName()   nogil except + # wrap-attach:FeatureFinderAlgorithmPicked
