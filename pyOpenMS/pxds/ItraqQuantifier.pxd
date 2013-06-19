from MSSpectrum cimport *
from FeatureMap cimport *
from ConsensusMap cimport *
from MSExperiment cimport *
from ChromatogramPeak cimport *
from Peak1D cimport *
from Param cimport *
from DefaultParamHandler cimport *
from ProgressLogger cimport *

from Types cimport *
from ItraqConstants cimport *
from ConsensusMap cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>" namespace "OpenMS":

    cdef cppclass ItraqQuantifier(ItraqConstants,DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler
        #    ItraqConstants

        ItraqQuantifier() nogil except +
        ItraqQuantifier(ItraqQuantifier) nogil except + #wrap-ignore
        ItraqQuantifier(Int itraq_type, Param param) nogil except +

        void run(ConsensusMap & map_in, ConsensusMap & map_out) nogil except +

        ItraqQuantifierStats getStats()


cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>" namespace "OpenMS::ItraqQuantifier":
    
    cdef cppclass ItraqQuantifierStats "OpenMS::ItraqQuantifier::ItraqQuantifierStats":
        ItraqQuantifierStats() nogil except +
        ItraqQuantifierStats(ItraqQuantifierStats) nogil except + #wrap-ignore
        Size channel_count
        Size iso_number_ms2_negative
        Size iso_number_reporter_negative
        Size iso_number_reporter_different
        DoubleReal iso_solution_different_intensity
        DoubleReal iso_total_intensity_negative
        Size number_ms2_total
        Size number_ms2_empty
        # TODO STL attribute
        libcpp_map[ size_t, size_t ] empty_channels

