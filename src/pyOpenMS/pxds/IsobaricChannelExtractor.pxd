from Types cimport *
from IsobaricQuantitationMethod cimport *
from ItraqFourPlexQuantitationMethod cimport *
from ItraqEightPlexQuantitationMethod cimport *
from TMTSixPlexQuantitationMethod cimport *
from TMTTenPlexQuantitationMethod cimport *
from DefaultParamHandler cimport *
from MSExperiment cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>" namespace "OpenMS":
    
    cdef cppclass IsobaricChannelExtractor(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler

        IsobaricChannelExtractor(IsobaricChannelExtractor &) nogil except +

        # IsobaricChannelExtractor(IsobaricQuantitationMethod *quant_method) nogil except + #wrap-ignore
        IsobaricChannelExtractor(ItraqEightPlexQuantitationMethod *quant_method) nogil except +
        IsobaricChannelExtractor(ItraqFourPlexQuantitationMethod *quant_method) nogil except +
        IsobaricChannelExtractor(TMTSixPlexQuantitationMethod *quant_method) nogil except +
        IsobaricChannelExtractor(TMTTenPlexQuantitationMethod *quant_method) nogil except +

        void extractChannels(MSExperiment & ms_exp_data, ConsensusMap & consensus_map) nogil except +

