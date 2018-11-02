from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from MSSpectrum cimport *
from MSExperiment cimport *
from Precursor cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PrecursorPurity.h>" namespace "OpenMS":

    cdef cppclass PrecursorPurity "OpenMS::PrecursorPurity":

        PrecursorPurity() nogil except +
        PrecursorPurity(PrecursorPurity) nogil except + # wrap-ignore

        libcpp_vector[PurityScores] computePrecursorPurities(MSExperiment spectra,
                                                             double precursor_mass_tolerance,
                                                             bool precursor_mass_tolerance_unit_ppm) nogil except +

        PurityScores computePrecursorPurity(MSSpectrum ms1,
                                            Precursor pre,
                                            double precursor_mass_tolerance,
                                            bool precursor_mass_tolerance_unit_ppm) nogil except +

    cdef cppclass PurityScores "OpenMS::PrecursorPurity::PurityScores":

        PurityScores(PurityScores) nogil except +
        PurityScores() nogil except +

        double total_intensity
        double target_intensity
        double residual_intensity
        double signal_proportion
        Size target_peak_count
        Size residual_peak_count
