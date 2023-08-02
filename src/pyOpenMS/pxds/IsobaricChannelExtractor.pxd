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

        IsobaricChannelExtractor(IsobaricChannelExtractor &) except + nogil 

        # IsobaricChannelExtractor(IsobaricQuantitationMethod *quant_method) except + nogil  #wrap-ignore
        IsobaricChannelExtractor(ItraqEightPlexQuantitationMethod *quant_method) except + nogil 
        IsobaricChannelExtractor(ItraqFourPlexQuantitationMethod *quant_method) except + nogil 
        IsobaricChannelExtractor(TMTSixPlexQuantitationMethod *quant_method) except + nogil 
        IsobaricChannelExtractor(TMTTenPlexQuantitationMethod *quant_method) except + nogil 

        void extractChannels(MSExperiment & ms_exp_data, ConsensusMap & consensus_map) except + nogil 

