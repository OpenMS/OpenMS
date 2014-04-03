from MSSpectrum cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>" namespace "OpenMS":

    cdef cppclass FeatureDeconvolution(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        FeatureDeconvolution()                  nogil except +
        FeatureDeconvolution(FeatureDeconvolution)   nogil except + #wrap-ignore

        void compute(FeatureMap[Feature] & input, FeatureMap[Feature] & output, ConsensusMap & cmap1, ConsensusMap & cmap2) nogil except +

cdef extern from "<OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>" namespace "OpenMS::FeatureDeconvolution":
    
    cdef enum CHARGEMODE "OpenMS::FeatureDeconvolution::CHARGEMODE":
        QFROMFEATURE
        QHEURISTIC
        QALL

