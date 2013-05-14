from ProgressLogger cimport *
from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from Types cimport *

from ConsensusMap cimport *
from MSExperiment cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/FILTERING/DATAREDUCTION/SILACAnalyzer.h>" namespace "OpenMS":

    cdef cppclass SILACAnalyzer(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        SILACAnalyzer()                 nogil except +
        SILACAnalyzer(SILACAnalyzer)    nogil except + # wrap-ignore

        # currently maps cannot be wrapped, do it manually... for wrap-ignore to work, it needs to be on one line!
        void initialize(String selected_labels_, UInt charge_min_, UInt charge_max_, Int missed_cleavages_, UInt isotopes_per_peptide_min_, UInt isotopes_per_peptide_max_, DoubleReal rt_threshold_, DoubleReal rt_min_, DoubleReal intensity_cutoff_, DoubleReal intensity_correlation_, DoubleReal model_deviation_, bool allow_missing_peaks_, libcpp_map[String, DoubleReal] label_identifiers) nogil except + # wrap-ignore

        void run_all(MSExperiment[Peak1D, ChromatogramPeak] & exp, ConsensusMap & cmap) nogil except +

        void writeConsensus(String & filename, ConsensusMap & cmap) nogil except +

