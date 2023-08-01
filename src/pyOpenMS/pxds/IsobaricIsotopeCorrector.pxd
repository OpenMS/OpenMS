from Types cimport *
from IsobaricQuantitationMethod cimport *
from ItraqFourPlexQuantitationMethod cimport *
from ItraqEightPlexQuantitationMethod cimport *
from TMTSixPlexQuantitationMethod cimport *
from TMTTenPlexQuantitationMethod cimport *
from ConsensusMap cimport *
from IsobaricQuantifierStatistics cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/IsobaricIsotopeCorrector.h>" namespace "OpenMS":

    cdef cppclass IsobaricIsotopeCorrector:

        IsobaricIsotopeCorrector() except + nogil  # compiler
        IsobaricIsotopeCorrector(IsobaricIsotopeCorrector &) except + nogil  # compiler

        # IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in,
        #                                                       ConsensusMap & consensus_map_out,
        #                                                       IsobaricQuantitationMethod * quant_method) except + nogil  # wrap-ignore

        IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in,
                                                               ConsensusMap & consensus_map_out,
                                                               ItraqEightPlexQuantitationMethod * quant_method) except + nogil 
        IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in,
                                                               ConsensusMap & consensus_map_out,
                                                               ItraqFourPlexQuantitationMethod * quant_method) except + nogil 
        IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in,
                                                               ConsensusMap & consensus_map_out,
                                                               TMTSixPlexQuantitationMethod * quant_method) except + nogil 
        IsobaricQuantifierStatistics correctIsotopicImpurities(ConsensusMap & consensus_map_in,
                                                               ConsensusMap & consensus_map_out,
                                                               TMTTenPlexQuantitationMethod * quant_method) except + nogil 

