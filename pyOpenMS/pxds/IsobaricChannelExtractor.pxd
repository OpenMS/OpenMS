from Types cimport *
from IsobaricQuantitationMethod cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>" namespace "OpenMS":
    
    cdef cppclass IsobaricChannelExtractor(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        IsobaricChannelExtractor(IsobaricChannelExtractor) nogil except +
        IsobaricChannelExtractor(IsobaricQuantitationMethod *quant_method) nogil except + #wrap-ignore
        void extractChannels(MSExperiment[ Peak1D, ChromatogramPeak ] & ms_exp_data, ConsensusMap & consensus_map) nogil except +

