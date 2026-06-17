from Types cimport *
from libcpp cimport bool
from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from String cimport *
from PeptideIdentification cimport *
from PeptideIdentificationList cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>" namespace "OpenMS":

    cdef cppclass OpenSearchModificationAnalysis:
        # wrap-doc:
        #  Utility class for analyzing modification patterns in open search results.
        #
        #  This class provides functionality to analyze delta mass patterns from open search
        #  peptide identifications, identify common modifications, and map them to known
        #  modifications from the ModificationsDB.
        #
        #  The class can generate two types of statistics tables:
        #  1. PTM Statistics Table - Shows known modifications mapped to delta masses
        #  2. Delta Mass Statistics Table - Shows raw delta mass distribution analysis
        #
        #  These tables can be used for modification discovery in open search workflows.
        #
        #  Usage:
        #
        #  .. code-block:: python
        #
        #    from pyopenms import *
        #    analyzer = OpenSearchModificationAnalysis()
        #
        #    # Load peptide identifications from open search
        #    peptides = []
        #    proteins = []
        #    IdXMLFile().load("open_search_results.idXML", proteins, peptides)
        #
        #    # Analyze modifications with statistics
        #    result = analyzer.analyzeModificationsWithStatistics(peptides, 5.0, True, False, "output.idXML")
        #
        #    # Access delta mass statistics
        #    dm_stats = result.delta_mass_stats
        #    print(f"Total PSMs: {dm_stats.total_psms}")
        #    print(f"Modified PSMs: {dm_stats.modified_psms}")
        #
        #    # Access PTM statistics
        #    ptm_stats = result.ptm_stats
        #    for ptm in ptm_stats.entries:
        #        print(f"{ptm.name}: {ptm.count} PSMs ({ptm.percentage}%)")

        OpenSearchModificationAnalysis() except + nogil
        OpenSearchModificationAnalysis(OpenSearchModificationAnalysis&) except + nogil

        # wrap-ignore (complex return type with custom comparator not supported by autowrap)
        # libcpp_pair[...] analyzeDeltaMassPatterns(...)

        # wrap-ignore (takes histogram types with custom comparator)
        # libcpp_vector[OpenSearchModificationAnalysis_ModificationSummary] mapDeltaMassesToModifications(...)

        libcpp_vector[OpenSearchModificationAnalysis_ModificationSummary] analyzeModifications(
            PeptideIdentificationList& peptide_ids,
            double precursor_mass_tolerance,
            bool precursor_mass_tolerance_unit_ppm,
            bool use_smoothing,
            String output_file) except + nogil
        # wrap-doc:
        #  Complete analysis workflow: analyze patterns and map to modifications.
        #
        #  :param peptide_ids: List of peptide identifications (modified in-place)
        #  :param precursor_mass_tolerance: Mass tolerance for mapping
        #  :param precursor_mass_tolerance_unit_ppm: Whether tolerance is in ppm
        #  :param use_smoothing: Apply smoothing to histogram
        #  :param output_file: Optional output file for results
        #  :returns: List of modification summaries

        OpenSearchModificationAnalysis_OpenSearchAnalysisResult analyzeModificationsWithStatistics(
            PeptideIdentificationList& peptide_ids,
            double precursor_mass_tolerance,
            bool precursor_mass_tolerance_unit_ppm,
            bool use_smoothing,
            String output_file) except + nogil
        # wrap-doc:
        #  Complete analysis returning structured statistics tables.
        #
        #  This is the main entry point for fragment index open search modification discovery.
        #  It performs a complete analysis workflow and returns structured tables.
        #
        #  :param peptide_ids: Peptide identifications (modified in-place with PTM annotations)
        #  :param precursor_mass_tolerance: Mass tolerance for mapping
        #  :param precursor_mass_tolerance_unit_ppm: Whether tolerance is in ppm
        #  :param use_smoothing: Apply Gaussian smoothing
        #  :param output_file: Optional output file path for TSV tables
        #  :returns: OpenSearchAnalysisResult with delta mass and PTM statistics

        # wrap-ignore (takes histogram types with custom comparator)
        # OpenSearchModificationAnalysis_DeltaMassStatistics generateDeltaMassStatistics(...)

        # wrap-ignore (takes histogram types with custom comparator)
        # OpenSearchModificationAnalysis_PTMStatistics generatePTMStatistics(...)

        libcpp_map[char, int] analyzeResidueFrequency(
            PeptideIdentificationList& peptide_ids,
            double delta_mass,
            double tolerance) except + nogil
        # wrap-doc:
        #  Analyze which amino acid residues are associated with a delta mass.
        #
        #  :param peptide_ids: Peptide identifications with DeltaMass meta values
        #  :param delta_mass: Target delta mass to analyze
        #  :param tolerance: Mass tolerance for matching
        #  :returns: Map from amino acid character to occurrence count

        void writeDeltaMassStatistics(
            OpenSearchModificationAnalysis_DeltaMassStatistics& stats,
            String output_file) except + nogil
        # wrap-doc:
        #  Write delta mass statistics to a TSV file.
        #
        #  :param stats: Delta mass statistics to write
        #  :param output_file: Output file path

        void writePTMStatistics(
            OpenSearchModificationAnalysis_PTMStatistics& stats,
            String output_file) except + nogil
        # wrap-doc:
        #  Write PTM statistics to a TSV file.
        #
        #  :param stats: PTM statistics to write
        #  :param output_file: Output file path


cdef extern from "<OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>" namespace "OpenMS::OpenSearchModificationAnalysis":

    cdef cppclass OpenSearchModificationAnalysis_ModificationPattern "OpenMS::OpenSearchModificationAnalysis::ModificationPattern":
        double count
        libcpp_vector[double] masses
        int num_charge_states

        OpenSearchModificationAnalysis_ModificationPattern() except + nogil
        OpenSearchModificationAnalysis_ModificationPattern(OpenSearchModificationAnalysis_ModificationPattern&) except + nogil

    cdef cppclass OpenSearchModificationAnalysis_ModificationSummary "OpenMS::OpenSearchModificationAnalysis::ModificationSummary":
        int count
        String name
        int num_charge_states
        libcpp_vector[double] masses

        OpenSearchModificationAnalysis_ModificationSummary() except + nogil
        OpenSearchModificationAnalysis_ModificationSummary(OpenSearchModificationAnalysis_ModificationSummary&) except + nogil

    cdef cppclass OpenSearchModificationAnalysis_DeltaMassEntry "OpenMS::OpenSearchModificationAnalysis::DeltaMassEntry":
        # wrap-doc:
        #  Statistics for a single delta mass bin in the histogram.
        double delta_mass
        int count
        int unique_peptides
        int num_charge_states
        double percentage
        String mapped_modification
        bool is_known_modification

        OpenSearchModificationAnalysis_DeltaMassEntry() except + nogil
        OpenSearchModificationAnalysis_DeltaMassEntry(OpenSearchModificationAnalysis_DeltaMassEntry&) except + nogil

    cdef cppclass OpenSearchModificationAnalysis_PTMEntry "OpenMS::OpenSearchModificationAnalysis::PTMEntry":
        # wrap-doc:
        #  Statistics for a mapped PTM.
        String name
        double theoretical_mass
        double observed_mass
        double mass_deviation
        int count
        int unique_peptides
        int num_charge_states
        double percentage
        # libcpp_map[char, int] not supported as struct member by autowrap;
        # use analyzeResidueFrequency() method for per-residue counts
        String target_residues

        OpenSearchModificationAnalysis_PTMEntry() except + nogil
        OpenSearchModificationAnalysis_PTMEntry(OpenSearchModificationAnalysis_PTMEntry&) except + nogil

    cdef cppclass OpenSearchModificationAnalysis_DeltaMassStatistics "OpenMS::OpenSearchModificationAnalysis::DeltaMassStatistics":
        # wrap-doc:
        #  Container for delta mass statistics table.
        libcpp_vector[OpenSearchModificationAnalysis_DeltaMassEntry] entries
        int total_psms
        int modified_psms
        int unmodified_psms
        double mean_delta_mass
        double median_delta_mass

        OpenSearchModificationAnalysis_DeltaMassStatistics() except + nogil
        OpenSearchModificationAnalysis_DeltaMassStatistics(OpenSearchModificationAnalysis_DeltaMassStatistics&) except + nogil

    cdef cppclass OpenSearchModificationAnalysis_PTMStatistics "OpenMS::OpenSearchModificationAnalysis::PTMStatistics":
        # wrap-doc:
        #  Container for PTM statistics table.
        libcpp_vector[OpenSearchModificationAnalysis_PTMEntry] entries
        int total_modified_psms
        int unknown_modification_psms
        int num_unique_modifications

        OpenSearchModificationAnalysis_PTMStatistics() except + nogil
        OpenSearchModificationAnalysis_PTMStatistics(OpenSearchModificationAnalysis_PTMStatistics&) except + nogil

    cdef cppclass OpenSearchModificationAnalysis_OpenSearchAnalysisResult "OpenMS::OpenSearchModificationAnalysis::OpenSearchAnalysisResult":
        # wrap-doc:
        #  Combined result of open search modification analysis.
        OpenSearchModificationAnalysis_DeltaMassStatistics delta_mass_stats
        OpenSearchModificationAnalysis_PTMStatistics ptm_stats
        libcpp_vector[OpenSearchModificationAnalysis_ModificationSummary] summaries

        OpenSearchModificationAnalysis_OpenSearchAnalysisResult() except + nogil
        OpenSearchModificationAnalysis_OpenSearchAnalysisResult(OpenSearchModificationAnalysis_OpenSearchAnalysisResult&) except + nogil
