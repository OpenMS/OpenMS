from Types cimport *
from IsobaricQuantitationMethod cimport *
from ItraqFourPlexQuantitationMethod cimport *
from ItraqEightPlexQuantitationMethod cimport *
from TMTSixPlexQuantitationMethod cimport *
from TMTTenPlexQuantitationMethod cimport *
from IsobaricQuantifierStatistics cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricNormalizer.h>" namespace "OpenMS":
    
    cdef cppclass IsobaricNormalizer "OpenMS::IsobaricNormalizer":
        IsobaricNormalizer(IsobaricNormalizer &) nogil except +

        IsobaricNormalizer(IsobaricQuantitationMethod *quant_method) nogil except + # wrap-ignore
        IsobaricNormalizer(ItraqFourPlexQuantitationMethod *quant_method) nogil except +
        IsobaricNormalizer(ItraqEightPlexQuantitationMethod *quant_method) nogil except +
        IsobaricNormalizer(TMTSixPlexQuantitationMethod *quant_method) nogil except +
        IsobaricNormalizer(TMTTenPlexQuantitationMethod *quant_method) nogil except +

        void normalize(ConsensusMap & consensus_map) nogil except +

