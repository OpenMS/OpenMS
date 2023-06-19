from Types cimport *
from DefaultParamHandler cimport *
from ConsensusMap cimport *
from IsobaricQuantitationMethod cimport *
from ItraqFourPlexQuantitationMethod cimport *
from ItraqEightPlexQuantitationMethod cimport *
from TMTSixPlexQuantitationMethod cimport *
from TMTTenPlexQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>" namespace "OpenMS":
    
    cdef cppclass IsobaricQuantifier(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        IsobaricQuantifier(IsobaricQuantifier &) nogil except +

        IsobaricQuantifier(IsobaricQuantitationMethod *quant_method) nogil except +  # wrap-ignore
        IsobaricQuantifier(ItraqFourPlexQuantitationMethod *quant_method) nogil except +
        IsobaricQuantifier(ItraqEightPlexQuantitationMethod *quant_method) nogil except +
        IsobaricQuantifier(TMTSixPlexQuantitationMethod *quant_method) nogil except +
        IsobaricQuantifier(TMTTenPlexQuantitationMethod *quant_method) nogil except +

        void quantify(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out) nogil except +

