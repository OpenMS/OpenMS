from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
# from libcpp.map cimport map as libcpp_map
from libcpp cimport bool
from MSSpectrum cimport *
from MSExperiment cimport *
from Precursor cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PrecursorPurity.h>" namespace "OpenMS":

    cdef cppclass PrecursorPurity "OpenMS::PrecursorPurity":
        # wrap-doc:
            #   Precursor purity or noise estimation
            #   -----
            #   This class computes metrics for precursor isolation window purity (or noise)
            #   The function extracts the peaks from an isolation window targeted for fragmentation
            #   and determines which peaks are isotopes of the target and which come from other sources
            #   The intensities of the assumed target peaks are summed up as the target intensity
            #   Using this information it calculates an intensity ratio for the relative intensity of the target
            #   compared to other sources
            #   These metrics are combined over the previous and the next MS1 spectrum

        PrecursorPurity() nogil except +
        PrecursorPurity(PrecursorPurity &) nogil except +

        # libcpp_map[String, PurityScores] computePrecursorPurities(MSExperiment spectra,
        #                                                     double precursor_mass_tolerance,
        #                                                     bool precursor_mass_tolerance_unit_ppm) nogil except +

        PurityScores computePrecursorPurity(MSSpectrum ms1,
                                            Precursor pre,
                                            double precursor_mass_tolerance,
                                            bool precursor_mass_tolerance_unit_ppm) nogil except +
            # wrap-doc:
            #   Compute precursor purity metrics for one MS2 precursor
            #   -----
            #   Note: This function is implemented in a general way and can also be used for e.g. MS3 precursor isolation windows in MS2 spectra
            #   -----
            #   :param ms1: The Spectrum containing the isolation window
            #   :param pre: The precursor containing the definition the isolation window
            #   :param precursor_mass_tolerance: The precursor tolerance. Is used for determining the targeted peak and deisotoping
            #   :param precursor_mass_tolerance_unit_ppm: The unit of the precursor tolerance

    cdef cppclass PurityScores "OpenMS::PrecursorPurity::PurityScores":

        PurityScores() nogil except +
        PurityScores(PurityScores &) nogil except +

        double total_intensity
        double target_intensity
        double signal_proportion
        Size target_peak_count
        Size residual_peak_count
