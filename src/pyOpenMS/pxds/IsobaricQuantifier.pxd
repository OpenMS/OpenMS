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
        # wrap-doc:
        #  Performs isotope correction and normalization of isobaric labeling data.
        #  Given extracted channel intensities, this class corrects for isotope
        #  impurities using a correction matrix and optionally normalizes the
        #  intensities for further downstream processing

        # wrap-inherits:
        #  DefaultParamHandler
        IsobaricQuantifier(IsobaricQuantifier &) except + nogil

        IsobaricQuantifier(IsobaricQuantitationMethod *quant_method) except + nogil   # wrap-ignore
        IsobaricQuantifier(ItraqFourPlexQuantitationMethod *quant_method) except + nogil  # wrap-doc:Constructor for iTRAQ 4-plex quantitation
        IsobaricQuantifier(ItraqEightPlexQuantitationMethod *quant_method) except + nogil  # wrap-doc:Constructor for iTRAQ 8-plex quantitation
        IsobaricQuantifier(TMTSixPlexQuantitationMethod *quant_method) except + nogil  # wrap-doc:Constructor for TMT 6-plex quantitation
        IsobaricQuantifier(TMTTenPlexQuantitationMethod *quant_method) except + nogil  # wrap-doc:Constructor for TMT 10-plex quantitation

        void quantify(ConsensusMap & consensus_map_in, ConsensusMap & consensus_map_out) except + nogil  # wrap-doc:Applies isotope correction and normalization to the raw isobaric channel intensities
