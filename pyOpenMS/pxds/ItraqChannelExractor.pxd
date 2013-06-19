from Types cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ItraqConstants cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>" namespace "OpenMS":

    cdef cppclass ItraqChannelExtractor(ItraqConstants,DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ItraqConstants

        ItraqChannelExtractor() nogil except +
        ItraqChannelExtractor(ItraqChannelExtractor) nogil except + #wrap-ignore
        ItraqChannelExtractor(Int itraq_type, Param param) nogil except +

        void run(MSExperiment[Peak1D, ChromatogramPeak] & ms_exp, ConsensusMap & map_out) nogil except +

