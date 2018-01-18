from Types cimport *
from libcpp cimport bool
from libcpp.pair cimport pair as libcpp_pair
from libcpp.vector cimport vector as libcpp_vector
from MSSpectrum cimport *
from Types cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/XQuestScores.h>" namespace "OpenMS":
    
    cdef cppclass XQuestScores "OpenMS::XQuestScores":
        XQuestScores(XQuestScores) nogil except + #wrap-ignore
        float preScore(Size matched_alpha, Size ions_alpha, Size matched_beta, Size ions_beta) nogil except +
        float preScore(Size matched_alpha, Size ions_alpha) nogil except +

        double matchOddsScore(MSSpectrum& theoretical_spec,
                              double fragment_mass_tolerance, 
                              bool fragment_mass_tolerance_unit_ppm, 
                              bool is_xlink_spectrum, 
                              Size n_charges) nogil except +


        double logOccupancyProb(MSSpectrum theoretical_spec,
                                Size matched_size,
                                double fragment_mass_tolerance,
                                bool fragment_mass_tolerance_unit_ppm) nogil except +
                              
                               

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

