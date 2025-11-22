from libcpp.vector cimport vector as libcpp_vector

from Types cimport *
from DefaultParamHandler cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from String cimport String

cdef extern from "<OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>" namespace "OpenMS":

    cdef cppclass Biosaur2Algorithm(DefaultParamHandler):
        cdef cppclass Hill:
            Hill() except + nogil
            libcpp_vector[Size] scan_indices
            libcpp_vector[Size] peak_indices
            libcpp_vector[double] mz_values
            libcpp_vector[double] intensities
            libcpp_vector[double] rt_values
            libcpp_vector[double] drift_times
            libcpp_vector[double] ion_mobilities
            double mz_median
            double rt_start
            double rt_end
            double rt_apex
            double intensity_apex
            double intensity_sum
            double drift_time_median
            double ion_mobility_median
            Size length
            Size hill_idx

        cdef cppclass IsotopeCandidate:
            IsotopeCandidate() except + nogil
            Size hill_idx
            Size isotope_number
            double mass_diff_ppm
            double cos_corr

        cdef cppclass PeptideFeature:
            PeptideFeature() except + nogil
            double mz
            double rt_start
            double rt_end
            double rt_apex
            double intensity_apex
            double intensity_sum
            int charge
            Size n_isotopes
            Size n_scans
            double mass_calib
            double drift_time
            double ion_mobility
            libcpp_vector[IsotopeCandidate] isotopes
            Size mono_hill_idx

        Biosaur2Algorithm() except + nogil

        void run(const MSExperiment& input,
                 FeatureMap& feature_map) except + nogil

        void run(const MSExperiment& input,
                 FeatureMap& feature_map,
                 libcpp_vector[Hill]& hills,
                 libcpp_vector[PeptideFeature]& features) except + nogil

        void writeTSV(const libcpp_vector[PeptideFeature]& features,
                      const String& filename) except + nogil

        void writeHills(const libcpp_vector[Hill]& hills,
                        const String& filename) except + nogil
