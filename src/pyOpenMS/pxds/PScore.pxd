from Types cimport *
from MSSpectrum cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/RNPXL/PScore.h>" namespace "OpenMS":
    
    cdef cppclass PScore "OpenMS::PScore":
        PScore(PScore) nogil except + #wrap-ignore

cdef extern from "<OpenMS/ANALYSIS/RNPXL/PScore.h>" namespace "OpenMS::PScore":

    libcpp_vector[ size_t ] calculateIntensityRankInMZWindow(libcpp_vector[ double ] & mz,
                                                             libcpp_vector[ double ] & intensities, 
                                                             double mz_window) nogil except + # wrap-attach:PScore

    libcpp_vector[ libcpp_vector[ size_t ] ] calculateRankMap(MSExperiment[Peak1D, ChromatogramPeak] & peak_map, 
                                                              double mz_window) nogil except + # wrap-attach:PScore

    # libcpp_map[ Size, MSSpectrum[Peak1D] ] calculatePeakLevelSpectra(MSSpectrum[Peak1D] & spec,
    #                                                                  libcpp_vector[ size_t ] & ranks,
    #                                                                  Size min_level,
    #                                                                  Size max_level) nogil except +

    double massCorrectionTerm(double mass) nogil except + # wrap-attach:PScore

    double cleavageCorrectionTerm(Size cleavages,
                                  bool consecutive_cleavage) nogil except + # wrap-attach:PScore


    double modificationCorrectionTerm(Size modifications) nogil except + # wrap-attach:PScore


    double computePScore(double fragment_mass_tolerance, 
                         bool fragment_mass_tolerance_unit_ppm,
                         libcpp_map[ size_t, MSSpectrum[Peak1D] ] & peak_level_spectra, 
                         libcpp_vector[ MSSpectrum[RichPeak1D] ] & theo_spectra, 
                         double mz_window) nogil except + # wrap-attach:PScore

    # overloaded static methods do not seem to work ...

    # double computePScore(double fragment_mass_tolerance, 
    #                      bool fragment_mass_tolerance_unit_ppm,
    #                      libcpp_map[ Size, MSSpectrum[Peak1D] ] & peak_level_spectra,
    #                      MSSpectrum[RichPeak1D] & theo_spectrum,
    #                      double mz_window) nogil except +

