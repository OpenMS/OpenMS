from Peak1D cimport *
from Feature cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from FeatureFinder cimport *
from DefaultParamHandler cimport *

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>" namespace "OpenMS":

    cdef cppclass FeatureFinderAlgorithmIsotopeWavelet[PeakT](DefaultParamHandler):

        # wrap-inherits:
        #    DefaultParamHandler
        #
        # wrap-instances:
        #    FeatureFinderAlgorithmIsotopeWavelet := FeatureFinderAlgorithmIsotopeWavelet[Peak1D]

        FeatureFinderAlgorithmIsotopeWavelet()      nogil except +

        void setData(MSExperiment[Peak1D, ChromatogramPeak] & input, FeatureMap& output, FeatureFinder & ff) nogil except +
        void run() nogil except +
        # MSSpectrum<PeakType> * createHRData(const UInt i)
        # static FeatureFinderAlgorithm<PeakType> * create()

cdef extern from "<OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmIsotopeWavelet.h>" namespace "OpenMS::FeatureFinderAlgorithmIsotopeWavelet<OpenMS::Peak1D>":

    String getProductName()   nogil except + # wrap-attach:FeatureFinderAlgorithmIsotopeWavelet
