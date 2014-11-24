from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmSH[PeakT](DefaultParamHandler):

        # wrap-inherits:
        #    DefaultParamHandler
        #
        # wrap-instances:
        #    FeatureFinderAlgorithmSH := FeatureFinderAlgorithmSH[Peak1D]

        FeatureFinderAlgorithmSH()      nogil except +

        void setData(MSExperiment[Peak1D, ChromatogramPeak] & input, FeatureMap & output, FeatureFinder & ff) nogil except +
        void run() nogil except +
        unsigned int getNativeScanId(String native_id) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>" namespace "OpenMS::FeatureFinderAlgorithmSH<OpenMS::Peak1D>":

    String getProductName()   nogil except + # wrap-attach:FeatureFinderAlgorithmSH
