from Types cimport *
from ConsensusMap cimport *
from IsobaricQuantitationMethod cimport *
from IsobaricQuantifierStatistics cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>" namespace "OpenMS":
    
    cdef cppclass IsobaricIsotopeCorrector "OpenMS::IsobaricIsotopeCorrector":
        IsobaricIsotopeCorrector(IsobaricIsotopeCorrector) nogil except +
        IsobaricIsotopeCorrector(IsobaricQuantitationMethod *quant_method) nogil except + # wrap-ignore
        IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out) nogil except +

