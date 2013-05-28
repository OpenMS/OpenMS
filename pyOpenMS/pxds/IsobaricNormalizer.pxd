from Types cimport *
from IsobaricQuantitationMethod cimport *
from IsobaricQuantifierStatistics cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricNormalizer.h>" namespace "OpenMS":
    
    cdef cppclass IsobaricNormalizer "OpenMS::IsobaricNormalizer":
        IsobaricNormalizer(IsobaricNormalizer) nogil except +
        IsobaricNormalizer(IsobaricQuantitationMethod *quant_method) nogil except + # wrap-ignore
        void normalize(ConsensusMap & consensus_map) nogil except +

