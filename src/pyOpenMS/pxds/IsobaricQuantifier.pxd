from Types cimport *
from DefaultParamHandler cimport *
from ConsensusMap cimport *
from IsobaricQuantitationMethod cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>" namespace "OpenMS":
    
    cdef cppclass IsobaricQuantifier(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        IsobaricQuantifier(IsobaricQuantifier) nogil except +
        IsobaricQuantifier(IsobaricQuantitationMethod *quant_method) nogil except +  # wrap-ignore
        void quantify(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out) nogil except +

