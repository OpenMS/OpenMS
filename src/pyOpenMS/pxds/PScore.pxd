from Types cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/PScore.h>" namespace "OpenMS":
    
    cdef cppclass PScore "OpenMS::PScore":
        PScore() nogil except + 
        PScore(PScore) nogil except + #wrap-ignore

        libcpp_vector[ size_t ] calculateIntensityRankInMZWindow(libcpp_vector[ double ] & mz, libcpp_vector[ double ] & intensities, double mz_window) nogil except +

        libcpp_vector[ libcpp_vector[ size_t ] ] calculateRankMap(MSExperiment & peak_map, double mz_window) nogil except +

        libcpp_map[ size_t, MSSpectrum ] calculatePeakLevelSpectra(MSSpectrum & spec, libcpp_vector[ size_t ] & ranks, Size min_level, Size max_level) nogil except +

        double computePScore(double fragment_mass_tolerance,
                             bool fragment_mass_tolerance_unit_ppm,
                             libcpp_map[ size_t, MSSpectrum ] & peak_level_spectra,
                             libcpp_vector[ MSSpectrum ] & theo_spectra,
                             double mz_window) nogil except +

        double computePScore(double fragment_mass_tolerance,
                             bool fragment_mass_tolerance_unit_ppm,
                             libcpp_map[ size_t, MSSpectrum ] & peak_level_spectra,
                             MSSpectrum & theo_spectrum,
                             double mz_window) nogil except +

## wrap static methods
# cdef extern from "<OpenMS/ANALYSIS/RNPXL/PScore.h>" namespace "OpenMS::PScore":
# 
#     double massCorrectionTerm(double mass) nogil except + #wrap-attach:PScore
#     double cleavageCorrectionTerm(Size cleavages, bool consecutive_cleavage) nogil except + #wrap-attach:PScore
#     double modificationCorrectionTerm(Size modifications) nogil except + #wrap-attach:PScore
# 
