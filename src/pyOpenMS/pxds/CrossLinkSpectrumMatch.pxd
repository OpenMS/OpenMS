from Types cimport *
from libcpp cimport bool
from PeptideHit cimport *
from MSExperiment cimport *
from ResidueModification cimport *
from FASTAFile cimport *
from ProteaseDigestion cimport *
from ProteinProteinCrossLink cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":
    
    cdef cppclass OPXL_CrossLinkSpectrumMatch "OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch":
        OPXL_CrossLinkSpectrumMatch(OPXL_CrossLinkSpectrumMatch) nogil except + #wrap-ignore
        ProteinProteinCrossLink cross_link
        Size scan_index_light
        Size scan_index_heavy
        double score
        Size rank
        double pre_score
        double percTIC
        double wTIC
        double int_sum
        double match_odds
        libcpp_vector[ double ] xcorrx
        double xcorrx_max
        libcpp_vector[ double ] xcorrc
        double xcorrc_max
        Size matched_common_alpha
        Size matched_common_beta
        Size matched_xlink_alpha
        Size matched_xlink_beta
        double HyperCommon
        double HyperXlink
        double HyperAlpha
        double HyperBeta
        double HyperBoth
        double PScoreCommon
        double PScoreXlink
        double PScoreAlpha
        double PScoreBeta
        double PScoreBoth
        libcpp_vector[ PeptideHit_PeakAnnotation ] frag_annotations
        Size peptide_id_index
        bool operator<(OPXL_CrossLinkSpectrumMatch & other) nogil except +
        bool operator==(OPXL_CrossLinkSpectrumMatch & other) nogil except +

