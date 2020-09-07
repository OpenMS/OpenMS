from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmSH(DefaultParamHandler):

        # wrap-inherits:
        #    DefaultParamHandler
        FeatureFinderAlgorithmSH()      nogil except +

        void setData(MSExperiment & input, FeatureMap & output, FeatureFinder & ff) nogil except +
        void run() nogil except +
        unsigned int getNativeScanId(String native_id) nogil except +

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>" namespace "OpenMS::FeatureFinderAlgorithmSH":

    String getProductName()   nogil except + # wrap-attach:FeatureFinderAlgorithmSH
