from OPXLDataStructs cimport *
from ProteinProteinCrossLink cimport ProteinProteinCrossLink
from PeptideHit cimport PeptideHit_PeakAnnotation

cdef extern from "<OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>" namespace "OpenMS::OPXLDataStructs":

    cdef cppclass CrossLinkSpectrumMatch "OpenMS::OPXLDataStructs::CrossLinkSpectrumMatch":

        CrossLinkSpectrumMatch(CrossLinkSpectrumMatch) nogil except +
        CrossLinkSpectrumMatch() nogil except +

        ProteinProteinCrossLink cross_link

        Size scan_index_light
        Size scan_index_heavy
        double score
        Size rank
        double xquest_score
        double pre_score
        double percTIC
        double wTIC
        double wTICold
        double int_sum
        double intsum_alpha
        double intsum_beta
        double total_current
        double precursor_error_ppm

        double match_odds
        double match_odds_alpha
        double match_odds_beta
        double log_occupancy
        double log_occupancy_alpha
        double log_occupancy_beta
        double xcorrx_max
        double xcorrc_max

        Size matched_linear_alpha
        Size matched_linear_beta
        Size matched_xlink_alpha
        Size matched_xlink_beta

        double num_iso_peaks_mean
        double num_iso_peaks_mean_linear_alpha
        double num_iso_peaks_mean_linear_beta
        double num_iso_peaks_mean_xlinks_alpha
        double num_iso_peaks_mean_xlinks_beta

        double ppm_error_abs_sum_linear_alpha
        double ppm_error_abs_sum_linear_beta
        double ppm_error_abs_sum_xlinks_alpha
        double ppm_error_abs_sum_xlinks_beta
        double ppm_error_abs_sum_linear
        double ppm_error_abs_sum_xlinks
        double ppm_error_abs_sum_alpha
        double ppm_error_abs_sum_beta
        double ppm_error_abs_sum

        int precursor_correction

        double precursor_total_intensity
        double precursor_target_intensity
        double precursor_signal_proportion
        Size precursor_target_peak_count
        Size precursor_residual_peak_count

        libcpp_vector[ PeptideHit_PeakAnnotation ] frag_annotations
        Size peptide_id_index
