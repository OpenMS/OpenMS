from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from MSSpectrum cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/XQuestScores.h>" namespace "OpenMS":

    cdef cppclass XQuestScores "OpenMS::XQuestScores":
        XQuestScores() nogil except +
        XQuestScores(XQuestScores &) nogil except +
        float preScore(Size matched_alpha, Size ions_alpha, Size matched_beta, Size ions_beta) nogil except +
        # wrap-doc:
                #   Compute a simple and fast to compute pre-score for a cross-link spectrum match
                #   -----
                #   :param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
                #   :param ions_alpha: Number of theoretical ions from the alpha peptide
                #   :param matched_beta: Number of experimental peaks matched to theoretical linear ions from the beta peptide
                #   :param ions_beta: Number of theoretical ions from the beta peptide

        float preScore(Size matched_alpha, Size ions_alpha) nogil except +
        # wrap-doc:
                #   Compute a simple and fast to compute pre-score for a mono-link spectrum match
                #   -----
                #   :param matched_alpha: Number of experimental peaks matched to theoretical linear ions from the alpha peptide
                #   :param ions_alpha: Number of theoretical ions from the alpha peptide

        double matchOddsScore(MSSpectrum& theoretical_spec,
                              double fragment_mass_tolerance,
                              bool fragment_mass_tolerance_unit_ppm,
                              bool is_xlink_spectrum,
                              Size n_charges) nogil except +
        # wrap-doc:
                #   Compute the match-odds score, a score based on the probability of getting the given number of matched peaks by chance
                #   -----
                #   :param theoretical_spec: Theoretical spectrum, sorted by position
                #   :param matched_size: Alignment between the theoretical and the experimental spectra
                #   :param fragment_mass_tolerance: Fragment mass tolerance of the alignment
                #   :param fragment_mass_tolerance_unit_ppm: Fragment mass tolerance unit of the alignment, true = ppm, false = Da
                #   :param is_xlink_spectrum: Type of cross-link, true = cross-link, false = mono-link
                #   :param n_charges: Number of considered charges in the theoretical spectrum

        double logOccupancyProb(MSSpectrum theoretical_spec,
                                Size matched_size,
                                double fragment_mass_tolerance,
                                bool fragment_mass_tolerance_unit_ppm) nogil except +
        # wrap-doc:
                #   Compute the logOccupancyProb score, similar to the match_odds, a score based on the probability of getting the given number of matched peaks by chance
                #   -----
                #   :param theoretical_spec: Theoretical spectrum, sorted by position
                #   :param matched_size: Number of matched peaks between experimental and theoretical spectra
                #   :param fragment_mass_tolerance: The tolerance of the alignment
                #   :param fragment_mass_tolerance_unit: The tolerance unit of the alignment, true = ppm, false = Da

        double weightedTICScoreXQuest(Size alpha_size, Size beta_size,
                                     double intsum_alpha, double intsum_beta,
                                     double total_current, bool type_is_cross_link) nogil except +

        double weightedTICScore(Size alpha_size, Size beta_size, double intsum_alpha,
                                double intsum_beta, double total_current,
                                bool type_is_cross_link) nogil except +

        double matchedCurrentChain(libcpp_vector[ libcpp_pair[ size_t, size_t ] ] & matched_spec_common,
                                  libcpp_vector[ libcpp_pair[ size_t, size_t ] ] & matched_spec_xlinks,
                                  MSSpectrum & spectrum_common_peaks,
                                  MSSpectrum & spectrum_xlink_peaks) nogil except +

        double totalMatchedCurrent(libcpp_vector[ libcpp_pair[ size_t, size_t ] ] & matched_spec_common_alpha,
                                   libcpp_vector[ libcpp_pair[ size_t, size_t ] ] & matched_spec_common_beta,
                                   libcpp_vector[ libcpp_pair[ size_t, size_t ] ] & matched_spec_xlinks_alpha,
                                   libcpp_vector[ libcpp_pair[ size_t, size_t ] ] & matched_spec_xlinks_beta,
                                   MSSpectrum & spectrum_common_peaks,
                                   MSSpectrum & spectrum_xlink_peaks) nogil except +

        libcpp_vector[ double ] xCorrelation(MSSpectrum & spec1,
                                             MSSpectrum & spec2,
                                             Int maxshift, double tolerance) nogil except +

        double xCorrelationPrescore(MSSpectrum & spec1, MSSpectrum & spec2,double tolerance) nogil except +
