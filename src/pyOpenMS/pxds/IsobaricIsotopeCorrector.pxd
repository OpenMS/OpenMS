from Types cimport *
from ConsensusMap cimport *
from IsobaricQuantitationMethod cimport *
from IsobaricQuantifierStatistics cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>" namespace "OpenMS":

    cdef cppclass IsobaricIsotopeCorrector:

        IsobaricIsotopeCorrector() nogil except +
        IsobaricIsotopeCorrector(IsobaricIsotopeCorrector) nogil except + # wrap-ignore

        IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out, IsobaricQuantitationMethod * quant_method) nogil except + # wrap-ignore

