from libcpp.vector cimport vector as libcpp_vector

from Types cimport *
from DefaultParamHandler cimport *
from FeatureMap cimport *
from MSExperiment cimport *
from String cimport String

cdef extern from "<OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>" namespace "OpenMS":

    cdef cppclass Biosaur2Algorithm(DefaultParamHandler):
        # wrap-doc:
        #  C++ implementation of the Biosaur2 feature detection workflow.
        cdef cppclass Hill:
            # wrap-doc:
            #  Container describing a single RT-contiguous hill.
            Hill() except + nogil  # wrap-doc:Construct an empty hill representation
            Hill(Hill&) except + nogil  # wrap-doc:Copy constructor
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
            # wrap-doc:
            #  Auxiliary structure linking a hill to a specific isotope rank.
            IsotopeCandidate() except + nogil  # wrap-doc:Construct an empty isotope candidate
            IsotopeCandidate(IsotopeCandidate&) except + nogil  # wrap-doc:Copy constructor
            Size hill_idx
            Size isotope_number
            double mass_diff_ppm
            double cos_corr

        cdef cppclass PeptideFeature:
            # wrap-doc:
            #  Aggregated peptide feature description produced by Biosaur2Algorithm.
            PeptideFeature() except + nogil  # wrap-doc:Construct an empty peptide feature container
            PeptideFeature(PeptideFeature&) except + nogil  # wrap-doc:Copy constructor
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

        Biosaur2Algorithm() except + nogil  # wrap-doc:Instantiate the Biosaur2 feature finding algorithm
        Biosaur2Algorithm(Biosaur2Algorithm&) except + nogil  # wrap-doc:Copy constructor

        void run(const MSExperiment& input,
                 FeatureMap& feature_map) except + nogil  # wrap-doc:Run the algorithm storing only the resulting features

        void run(const MSExperiment& input,
                 FeatureMap& feature_map,
                 libcpp_vector[Hill]& hills,
                 libcpp_vector[PeptideFeature]& features) except + nogil  # wrap-doc:Run the algorithm and capture intermediate hills as well as peptide feature records

        void writeTSV(const libcpp_vector[PeptideFeature]& features,
                      const String& filename) except + nogil  # wrap-doc:Write detected peptide features to a Biosaur2-compatible TSV file

        void writeHills(const libcpp_vector[Hill]& hills,
                        const String& filename) except + nogil  # wrap-doc:Export the detected hills into a TSV diagnostic table
