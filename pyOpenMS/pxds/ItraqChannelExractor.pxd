from MSSpectrum cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>" namespace "OpenMS":

    cdef cppclass ItraqChannelExtractor(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        ItraqChannelExtractor()                  nogil except +
        ItraqChannelExtractor(ItraqChannelExtractor)   nogil except + #wrap-ignore

        void run(MSExperiment[Peak1D, ChromatogramPeak] & ms_exp, ConsensusMap & map_out) nogil except +

