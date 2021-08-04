from Types cimport *
from String cimport *
from ProgressLogger cimport *
from TargetedExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/OPENSWATH/MRMDecoy.h>" namespace "OpenMS":

    cdef cppclass MRMDecoy(ProgressLogger):
        # wrap-inherits:
        #    ProgressLogger

        MRMDecoy() nogil except +
        MRMDecoy(MRMDecoy &) nogil except + # compiler

        void generateDecoys(TargetedExperiment& exp,
                            TargetedExperiment& dec,
                            String method,
                            double aim_decoy_fraction,
                            bool switchKR,
                            String decoy_tag,
                            int max_attempts,
                            double identity_threshold,
                            double precursor_mz_shift,
                            double product_mz_shift,
                            double product_mz_threshold,
                            libcpp_vector[String] fragment_types,
                            libcpp_vector[size_t] fragment_charges,
                            bool enable_specific_losses,
                            bool enable_unspecific_losses,
                            int round_decPow) nogil except +
            # wrap-doc:
                #   Generate decoys from a TargetedExperiment
                #   -----
                #   Will generate decoy peptides for each target peptide provided in exp and
                #   write them into the decoy experiment
                #   -----
                #   Valid methods: shuffle, reverse, pseudo-reverse
                #   -----
                #   If theoretical is true, the target transitions will be returned but their
                #   masses will be adjusted to match the theoretical value of the fragment ion
                #   that is the most likely explanation for the product
                #   -----
                #   `mz_threshold` is used for the matching of theoretical ion series to the observed one
                #   -----
                #   To generate decoys with different precursor mass, use the "switchKR" flag
                #   which switches terminal K/R (switches K to R and R to K). This generates
                #   different precursor m/z and ensures that the y ion series has a different
                #   mass. For a description of the procedure, see (supplemental material)
                #   -----
                #   Bruderer et al. Mol Cell Proteomics. 2017. 10.1074/mcp.RA117.000314.


        libcpp_vector[size_t] findFixedResidues(const String & sequence,
                                                bool keepN,
                                                bool keepC,
                                                const String & keep_const_pattern) nogil except +
            # wrap-doc:
                #   Find all residues in a sequence that should not be reversed / shuffled
                #   -----
                #   :param sequence: The amino acid sequence
                #   :param keepN: Whether to keep N terminus constant
                #   :param keepC: Whether to keep C terminus constant
                #   :param keep_const_pattern: A string containing the AA to not change (e.g. 'KRP')

