from libcpp.vector cimport vector as libcpp_vector
from libcpp.map cimport map as libcpp_map
from libcpp cimport bool

from Types cimport *
from String cimport *
from PeptideIdentification cimport *
from PeptideIdentificationList cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>" namespace "OpenMS::OpenSearchModificationAnalysis":

    cdef cppclass OpenSearchModificationAnalysis_ModificationPattern "OpenMS::OpenSearchModificationAnalysis::ModificationPattern":
        # wrap-attach:
        #   OpenSearchModificationAnalysis
        # wrap-doc:
        #  Stores details of a modification pattern found in the data.

        OpenSearchModificationAnalysis_ModificationPattern() except + nogil  # compiler
        OpenSearchModificationAnalysis_ModificationPattern(OpenSearchModificationAnalysis_ModificationPattern &) except + nogil  # compiler

        double count
        libcpp_vector[double] masses
        int num_charge_states

    cdef cppclass OpenSearchModificationAnalysis_ModificationSummary "OpenMS::OpenSearchModificationAnalysis::ModificationSummary":
        # wrap-attach:
        #   OpenSearchModificationAnalysis
        # wrap-doc:
        #  Summary of a modification found in the data.
        #
        #  Attributes:
        #    count: Number of occurrences
        #    name: Modification name
        #    num_charge_states: Number of charge states observed
        #    masses: Associated masses

        OpenSearchModificationAnalysis_ModificationSummary() except + nogil  # compiler
        OpenSearchModificationAnalysis_ModificationSummary(OpenSearchModificationAnalysis_ModificationSummary &) except + nogil  # compiler

        int count
        String name
        int num_charge_states
        libcpp_vector[double] masses

cdef extern from "<OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>" namespace "OpenMS":

    cdef cppclass OpenSearchModificationAnalysis:
        # wrap-doc:
        #  Utility class for analyzing modification patterns in open search results.
        #
        #  Analyzes delta mass patterns from open search peptide identifications,
        #  identifies common modifications, and maps them to known modifications
        #  from the ModificationsDB.
        #
        #  Usage:
        #
        #  .. code-block:: python
        #
        #    from pyopenms import OpenSearchModificationAnalysis
        #
        #    analyzer = OpenSearchModificationAnalysis()
        #    summaries = analyzer.analyzeModifications(
        #      peptide_ids,
        #      precursor_mass_tolerance=0.01,
        #      precursor_mass_tolerance_unit_ppm=False,
        #      use_smoothing=True
        #    )
        #    for mod in summaries:
        #      print(mod.name, mod.count, mod.masses)

        OpenSearchModificationAnalysis() except + nogil  # wrap-doc:Default constructor
        OpenSearchModificationAnalysis(OpenSearchModificationAnalysis &) except + nogil  # compiler

        libcpp_vector[OpenSearchModificationAnalysis_ModificationSummary] analyzeModifications(
            PeptideIdentificationList & peptide_ids,
            double precursor_mass_tolerance,
            bool precursor_mass_tolerance_unit_ppm,
            bool use_smoothing,
            const String & output_file) except + nogil
            # wrap-doc:
            #  Complete analysis workflow: analyze delta mass patterns and map to modifications.
            #  :param peptide_ids: Peptide identifications with DeltaMass metavalues (modified in-place with PTM annotations)
            #  :param precursor_mass_tolerance: Mass tolerance for mapping delta masses to known modifications
            #  :param precursor_mass_tolerance_unit_ppm: True if tolerance is in ppm, False for Da
            #  :param use_smoothing: Whether to smooth the delta mass histogram before peak finding
            #  :param output_file: Optional file path for writing a modification summary table (empty string to skip)
            #  :returns: List of ModificationSummary objects sorted by relevance
