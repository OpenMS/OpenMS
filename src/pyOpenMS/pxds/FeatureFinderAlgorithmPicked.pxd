from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmPicked[PeakT](DefaultParamHandler):

        # wrap-inherits:
        #    DefaultParamHandler
        #
        # wrap-instances:
        #    FeatureFinderAlgorithmPicked := FeatureFinderAlgorithmPicked[Peak1D]

        FeatureFinderAlgorithmPicked()      nogil except +

        void setData(MSExperiment[Peak1D, ChromatogramPeak] & input, FeatureMap & output, FeatureFinder & ff) nogil except +
        void run() nogil except +


        void setSeeds(FeatureMap& seeds) nogil except +

        # static FeatureFinderAlgorithm<PeakType>* create()

#
# static methods are wrapped like this:
#

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPicked.h>" namespace "OpenMS::FeatureFinderAlgorithmPicked<OpenMS::Peak1D>":

    String getProductName()   nogil except + # wrap-attach:FeatureFinderAlgorithmPicked
