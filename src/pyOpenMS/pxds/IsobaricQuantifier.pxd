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
        IsobaricQuantifier(IsobaricQuantifier &) except + nogil 

        IsobaricQuantifier(IsobaricQuantitationMethod *quant_method) except + nogil   # wrap-ignore
        IsobaricQuantifier(ItraqFourPlexQuantitationMethod *quant_method) except + nogil 
        IsobaricQuantifier(ItraqEightPlexQuantitationMethod *quant_method) except + nogil 
        IsobaricQuantifier(TMTSixPlexQuantitationMethod *quant_method) except + nogil 
        IsobaricQuantifier(TMTTenPlexQuantitationMethod *quant_method) except + nogil 

        void quantify(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out) except + nogil 

