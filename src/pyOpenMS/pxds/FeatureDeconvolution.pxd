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
        #   DefaultParamHandler

        FeatureDeconvolution() except + nogil 
        FeatureDeconvolution(FeatureDeconvolution&) except + nogil  

        void compute(FeatureMap & input, FeatureMap & output, ConsensusMap & cmap1, ConsensusMap & cmap2) except + nogil 

cdef extern from "<OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>" namespace "OpenMS::FeatureDeconvolution":

    cdef enum class CHARGEMODE "OpenMS::FeatureDeconvolution::CHARGEMODE":
        #wrap-attach:
        #   FeatureDeconvolution
        QFROMFEATURE
        QHEURISTIC
        QALL

