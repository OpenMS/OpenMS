from Types cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/PScore.h>" namespace "OpenMS":
    
    cdef cppclass PScore "OpenMS::PScore":
        PScore() nogil except + 
        PScore(PScore &) nogil except + # compiler

        libcpp_vector[ size_t ] calculateIntensityRankInMZWindow(libcpp_vector[ double ] & mz, libcpp_vector[ double ] & intensities, double mz_window) nogil except +
            # wrap-doc:
                #   Calculate local (windowed) peak ranks
                #   -----
                #   The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity 
                #   The result can be used to efficiently filter spectra for top 1..n peaks in mass windows
                #   -----
                #   :param mz: The m/z positions of the peaks
                #   :param intensities: The intensities of the peaks
                #   :param mz_window: The window in Thomson centered at each peak 

        libcpp_vector[ libcpp_vector[ size_t ] ] calculateRankMap(MSExperiment & peak_map, double mz_window) nogil except +
            # wrap-doc:
                #   Precalculated, windowed peak ranks for a whole experiment
                #   -----
                #   The peak rank is defined as the number of neighboring peaks in +/- (mz_window/2) that have higher intensity 
                #   -----
                #   :param peak_map: Fragment spectra used for rank calculation. Typically a peak map after removal of all MS1 spectra
                #   :param mz_window: Window in Thomson centered at each peak

        libcpp_map[ size_t, MSSpectrum ] calculatePeakLevelSpectra(MSSpectrum & spec, libcpp_vector[ size_t ] & ranks, Size min_level, Size max_level) nogil except +
            # wrap-doc:
                #   Calculates spectra for peak level between min_level to max_level and stores them in the map
                #   -----
                #   A spectrum of peak level n retains the (n+1) top intensity peaks in a sliding mz_window centered at each peak

        double computePScore(double fragment_mass_tolerance,
                             bool fragment_mass_tolerance_unit_ppm,
                             libcpp_map[ size_t, MSSpectrum ] & peak_level_spectra,
                             libcpp_vector[ MSSpectrum ] & theo_spectra,
                             double mz_window) nogil except +
            # wrap-doc:
                #   Computes the PScore for a vector of theoretical spectra
                #   -----
                #   Similar to Andromeda, a vector of theoretical spectra can be provided that e.g. contain loss spectra or higher charge spectra depending on the sequence.
                #   The best score obtained by scoring all those theoretical spectra against the experimental ones is returned
                #   -----
                #   :param fragment_mass_tolerance: Mass tolerance for matching peaks
                #   :param fragment_mass_tolerance_unit_ppm: Whether Thomson or ppm is used
                #   :param peak_level_spectra: Spectra for different peak levels (=filtered by maximum rank).
                #   :param theo_spectra: Theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
                #   :param mz_window: Window in Thomson centered at each peak

        double computePScore(double fragment_mass_tolerance,
                             bool fragment_mass_tolerance_unit_ppm,
                             libcpp_map[ size_t, MSSpectrum ] & peak_level_spectra,
                             MSSpectrum & theo_spectrum,
                             double mz_window) nogil except +
            # wrap-doc:
                #   Computes the PScore for a single theoretical spectrum
                #   -----
                #   :param fragment_mass_tolerance: Mass tolerance for matching peaks
                #   :param fragment_mass_tolerance_unit_ppm: Whether Thomson or ppm is used
                #   :param peak_level_spectra: Spectra for different peak levels (=filtered by maximum rank)
                #   :param theo_spectra: Theoretical spectra as obtained e.g. from TheoreticalSpectrumGenerator
                #   :param mz_window: Window in Thomson centered at each peak

## wrap static methods
# cdef extern from "<OpenMS/ANALYSIS/RNPXL/PScore.h>" namespace "OpenMS::PScore":
# 
#     double massCorrectionTerm(double mass) nogil except + #wrap-attach:PScore
#     double cleavageCorrectionTerm(Size cleavages, bool consecutive_cleavage) nogil except + #wrap-attach:PScore
#     double modificationCorrectionTerm(Size modifications) nogil except + #wrap-attach:PScore
# 
