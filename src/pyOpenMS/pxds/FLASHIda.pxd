from libcpp.vector cimport vector as libcpp_vector
from libcpp.string cimport string as libcpp_string
from libcpp cimport bool

from String cimport *
from FASTAFile cimport *
from Param cimport *
from FLASHHelperClasses cimport *

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>" namespace "OpenMS":

    cdef cppclass FLASHIda:
        # wrap-doc:
        #   FLASHIda class for real-time deconvolution and proteoform identification.
        #
        #   This class contains functions to perform deconvolution for spectra received
        #   from mass spectrometers, as well as proteoform identification using the
        #   FLASHTnT workflow.

        # Constructor taking a string argument
        FLASHIda(char* arg) except + nogil

        # Copy constructor
        FLASHIda(FLASHIda &) except + nogil

        # Python-friendly method for sequence tag generation and matching
        int getSequenceTagsAndMatchesPy(libcpp_vector[double] & mzs,
                                        libcpp_vector[double] & ints,
                                        double rt,
                                        int ms_level,
                                        libcpp_vector[FASTAEntry] & fasta_entries,
                                        Param & tagger_param,
                                        libcpp_vector[Tag] & tags,
                                        libcpp_vector[FLASHIdaTagMatch] & matches,
                                        double ppm_tolerance,
                                        double max_flanking_mass_diff) except + nogil
            # wrap-doc:
            #   Deconvolute a spectrum and find sequence tags with database matches.
            #
            #   This method performs spectrum deconvolution using the FLASHDeconv algorithm,
            #   generates sequence tags using FLASHTagger, and matches them against a protein database.
            #
            #   :param mzs: m/z values of the input spectrum
            #   :param ints: intensities of the input spectrum
            #   :param rt: retention time in seconds
            #   :param ms_level: MS level of the spectrum
            #   :param fasta_entries: protein database entries to match against
            #   :param tagger_param: parameters for the FLASHTagger algorithm
            #   :param tags: output vector of detected sequence tags
            #   :param matches: output vector of tag matches to database entries
            #   :param ppm_tolerance: mass tolerance in ppm for tag matching
            #   :param max_flanking_mass_diff: maximum allowed flanking mass difference
            #   :returns: number of tags found

        # Python-friendly proteoform identification
        int identifyProteoformPy(libcpp_vector[double] & mzs,
                                 libcpp_vector[double] & ints,
                                 double rt,
                                 String & protein_sequence,
                                 double ppm_tolerance,
                                 libcpp_vector[String] & ion_types,
                                 double ptm_mass_threshold,
                                 libcpp_vector[int] & matched_fragment_indices,
                                 libcpp_vector[int] & ptm_start_positions,
                                 libcpp_vector[int] & ptm_end_positions,
                                 libcpp_vector[double] & ptm_masses) except + nogil
            # wrap-doc:
            #   Identify proteoform from MS2 spectrum against a single protein sequence.
            #
            #   Implements the core FLASHTnT identification workflow:
            #   1. Deconvolves the MS2 spectrum to get monoisotopic masses
            #   2. Generates sequence tags using FLASHTagger
            #   3. Runs FLASHExtender for PTM detection via DAG path-finding
            #
            #   :param mzs: m/z values of the input MS2 spectrum
            #   :param ints: intensities of the input MS2 spectrum
            #   :param rt: retention time in seconds
            #   :param protein_sequence: the protein sequence to match against
            #   :param ppm_tolerance: mass tolerance in ppm
            #   :param ion_types: ion types to consider (e.g., ["b", "y"])
            #   :param ptm_mass_threshold: (unused - PTM detection handled by FLASHExtender)
            #   :param matched_fragment_indices: output indices of matched fragment ions (1-based)
            #   :param ptm_start_positions: output start positions of PTM localization ranges (1-based)
            #   :param ptm_end_positions: output end positions of PTM localization ranges (1-based)
            #   :param ptm_masses: output mass shifts at each PTM position
            #   :returns: number of matched fragment ions

        # Python-friendly extended proteoform identification
        int identifyProteoformExtendedPy(libcpp_vector[double] & mzs,
                                         libcpp_vector[double] & ints,
                                         double rt,
                                         String & protein_sequence,
                                         double ppm_tolerance,
                                         libcpp_vector[String] & ion_types,
                                         int max_ptm_count,
                                         double max_ptm_mass,
                                         libcpp_vector[int] & matched_peak_indices,
                                         libcpp_vector[double] & matched_theoretical_masses,
                                         libcpp_vector[bool] & matched_ion_types,
                                         libcpp_vector[int] & ptm_start_positions,
                                         libcpp_vector[int] & ptm_end_positions,
                                         libcpp_vector[double] & ptm_masses,
                                         double & coverage,
                                         double & total_score) except + nogil
            # wrap-doc:
            #   Extended proteoform identification with detailed output.
            #
            #   Extended version providing more detailed match information including
            #   peak indices, theoretical masses, ion types, and PTM localization ranges.
            #
            #   :param mzs: m/z values of the input MS2 spectrum
            #   :param ints: intensities of the input MS2 spectrum
            #   :param rt: retention time in seconds
            #   :param protein_sequence: the protein sequence to match against
            #   :param ppm_tolerance: mass tolerance in ppm
            #   :param ion_types: ion types to consider
            #   :param max_ptm_count: maximum number of PTMs to consider
            #   :param max_ptm_mass: maximum mass shift for a single PTM
            #   :param matched_peak_indices: output indices of matched peaks
            #   :param matched_theoretical_masses: output theoretical masses matched
            #   :param matched_ion_types: output ion types (True=prefix, False=suffix)
            #   :param ptm_start_positions: output start of region for each PTM
            #   :param ptm_end_positions: output end of region for each PTM
            #   :param ptm_masses: output mass shift for each PTM
            #   :param coverage: output sequence coverage (0.0-1.0)
            #   :param total_score: output total identification score
            #   :returns: number of matched ions

        # Python-friendly identification from pre-deconvolved masses
        int identifyProteoformFromMassesPy(libcpp_vector[double] & observed_masses,
                                           libcpp_vector[double] & mass_scores,
                                           String & protein_sequence,
                                           double ppm_tolerance,
                                           libcpp_vector[String] & ion_types,
                                           double ptm_mass_threshold,
                                           libcpp_vector[int] & matched_fragment_indices,
                                           libcpp_vector[bool] & matched_ion_types,
                                           libcpp_vector[double] & matched_observed_masses,
                                           libcpp_vector[double] & matched_theoretical_masses,
                                           libcpp_vector[double] & matched_ppm_errors,
                                           libcpp_vector[int] & ptm_start_positions,
                                           libcpp_vector[int] & ptm_end_positions,
                                           libcpp_vector[double] & ptm_masses) except + nogil
            # wrap-doc:
            #   Identify proteoform from pre-deconvolved masses.
            #
            #   This function performs fragment ion matching without deconvolution,
            #   allowing direct testing with known masses. It takes deconvolved
            #   monoisotopic masses directly and matches them against theoretical
            #   fragment masses calculated from the protein sequence.
            #
            #   :param observed_masses: deconvolved monoisotopic masses
            #   :param mass_scores: quality scores for each mass (use 1.0 if unknown)
            #   :param protein_sequence: the protein sequence to match against
            #   :param ppm_tolerance: mass tolerance in ppm
            #   :param ion_types: ion types to consider (e.g., ["b", "y"])
            #   :param ptm_mass_threshold: (unused - PTM detection handled by FLASHExtender)
            #   :param matched_fragment_indices: output indices of matched fragments (1-based)
            #   :param matched_ion_types: output True=prefix/b-ion, False=suffix/y-ion
            #   :param matched_observed_masses: output observed mass for each match
            #   :param matched_theoretical_masses: output theoretical mass for each match
            #   :param matched_ppm_errors: output ppm error for each match
            #   :param ptm_start_positions: output start positions of PTM localization ranges (1-based)
            #   :param ptm_end_positions: output end positions of PTM localization ranges (1-based)
            #   :param ptm_masses: output mass shifts at PTM positions
            #   :returns: number of matched fragment ions

        # Static method for calculating theoretical fragment masses
        void calculateTheoreticalFragmentMassesPy(String & sequence,
                                                  libcpp_vector[String] & ion_types,
                                                  libcpp_vector[double] & prefix_masses,
                                                  libcpp_vector[double] & suffix_masses) except + nogil  # wrap-attach:FLASHIda
            # wrap-doc:
            #   Calculate theoretical fragment masses for a protein sequence.
            #
            #   :param sequence: the protein sequence
            #   :param ion_types: ion types to consider (e.g., ["b", "y"])
            #   :param prefix_masses: output cumulative masses from N-terminus
            #   :param suffix_masses: output cumulative masses from C-terminus

cdef extern from "<OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>" namespace "OpenMS::FLASHIda":

    cdef cppclass FLASHIdaTagMatch "OpenMS::FLASHIda::TagMatch":
        # wrap-doc:
        #   Structure representing a tag match to a protein database entry.

        FLASHIdaTagMatch() except + nogil
        FLASHIdaTagMatch(FLASHIdaTagMatch &) except + nogil

        String tag_sequence
        double n_term_mass
        double c_term_mass
        double tag_score
        int protein_index
        String protein_accession
        int match_position
        double flanking_mass_diff
