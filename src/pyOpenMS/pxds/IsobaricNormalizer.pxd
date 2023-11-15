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
        IsobaricNormalizer(IsobaricNormalizer &) except + nogil 

        IsobaricNormalizer(IsobaricQuantitationMethod *quant_method) except + nogil  # wrap-ignore
        IsobaricNormalizer(ItraqFourPlexQuantitationMethod *quant_method) except + nogil 
        IsobaricNormalizer(ItraqEightPlexQuantitationMethod *quant_method) except + nogil 
        IsobaricNormalizer(TMTSixPlexQuantitationMethod *quant_method) except + nogil 
        IsobaricNormalizer(TMTTenPlexQuantitationMethod *quant_method) except + nogil 

        void normalize(ConsensusMap & consensus_map) except + nogil 

