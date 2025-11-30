from libcpp.vector cimport vector as libcpp_vector
from libcpp cimport bool
from DefaultParamHandler cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from String cimport *
from ProgressLogger cimport *
from OpenSearchModificationAnalysis cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h>" namespace "OpenMS":

    cdef cppclass PeptideSearchEngineFIAlgorithm(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #  DefaultParamHandler
        #  ProgressLogger
        # wrap-doc:
        #   Fragment-index-based peptide database search algorithm (experimental).
        #
        #   Provides a self-contained search engine that matches MS/MS spectra against a protein
        #   database using an FI (Fragment Index). Typical usage:
        #   - Configure parameters via DefaultParamHandler (mass tolerances, enzyme, charges, etc.)
        #   - Call search() with an input mzML file and a FASTA database to populate identification
        #     outputs (ProteinIdentification and PeptideIdentificationList)
        #   - For modification discovery in open search mode, use searchWithModificationAnalysis()
        #     which returns structured PTM and delta mass statistics tables
        #   - Intended for educational/prototyping use and to demonstrate FI-backed searching
        #
        #   Usage for open search modification discovery:
        #
        #   .. code-block:: python
        #
        #     from pyopenms import *
        #     algo = PeptideSearchEngineFIAlgorithm()
        #     p = algo.getParameters()
        #     p.setValue("precursor:mass_tolerance", 500.0)  # Open search
        #     p.setValue("precursor:mass_tolerance_unit", "Da")
        #     algo.setParameters(p)
        #
        #     result = algo.searchWithModificationAnalysis("spectra.mzML", "database.fasta", "output")
        #     if result.exit_code == PeptideSearchEngineFIAlgorithm.EXECUTION_OK and result.is_open_search:
        #         # Access PTM statistics
        #         for ptm in result.modification_analysis.ptm_stats.entries:
        #             print(f"{ptm.name}: {ptm.count} PSMs")

        PeptideSearchEngineFIAlgorithm() except + nogil  # compiler
        PeptideSearchEngineFIAlgorithm(PeptideSearchEngineFIAlgorithm &) except + nogil  # compiler

        PeptideSearchEngineFIAlgorithm_ExitCodes search(const String& in_mzML,
          const String& in_db,
          libcpp_vector[ ProteinIdentification ] & prot_ids,
          PeptideIdentificationList & pep_ids) except + nogil
        # wrap-doc:
        #   Search spectra in an mzML file against a protein database using an FI-backed workflow.
        #
        #   Populates protein and peptide identifications, including search meta data, PSM hits,
        #   and search engine annotations. Parameters are taken from this instance.
        #
        #   :param in_mzML: Input path to the mzML file containing MS/MS spectra to search
        #   :param in_db: Input path to the protein sequence database in FASTA format
        #   :param prot_ids: Output container receiving search meta data and protein-level information
        #   :param pep_ids: Output container receiving spectrum-level peptide identifications (PSMs)
        #   :returns: ExitCodes indicating success (EXECUTION_OK) or the encountered error condition

        PeptideSearchEngineFIAlgorithm_SearchResult searchWithModificationAnalysis(
          String in_mzML,
          String in_db,
          String output_base_name) except + nogil
        # wrap-doc:
        #   Search with comprehensive results including modification analysis tables.
        #
        #   This method performs a peptide database search and additionally returns
        #   structured modification analysis results for open search mode. This is the
        #   recommended method for modification discovery workflows.
        #
        #   The method automatically:
        #   - Detects open search mode based on precursor tolerance
        #   - Computes delta mass statistics
        #   - Maps delta masses to known modifications
        #   - Generates PTM statistics with residue localization
        #   - Writes TSV output files if output_base_name is provided
        #
        #   :param in_mzML: Input path to the mzML file
        #   :param in_db: Input path to the FASTA database
        #   :param output_base_name: Optional base name for output files (TSV tables)
        #   :returns: SearchResult containing identifications and modification analysis

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h>" namespace "OpenMS::PeptideSearchEngineFIAlgorithm":
    cdef enum PeptideSearchEngineFIAlgorithm_ExitCodes "OpenMS::PeptideSearchEngineFIAlgorithm::ExitCodes":
        #wrap-attach:
        #   PeptideSearchEngineFIAlgorithm
        EXECUTION_OK
        INPUT_FILE_EMPTY
        UNEXPECTED_RESULT
        UNKNOWN_ERROR
        ILLEGAL_PARAMETERS

    cdef cppclass PeptideSearchEngineFIAlgorithm_SearchResult "OpenMS::PeptideSearchEngineFIAlgorithm::SearchResult":
        # wrap-attach:
        #   PeptideSearchEngineFIAlgorithm
        # wrap-doc:
        #   Comprehensive search result including modification analysis.
        #
        #   This structure contains all outputs from an open search including:
        #   - Standard protein and peptide identifications
        #   - Delta mass statistics table (histogram of mass shifts)
        #   - PTM statistics table (mapped modifications with residue analysis)
        PeptideSearchEngineFIAlgorithm_ExitCodes exit_code
        libcpp_vector[ProteinIdentification] protein_ids
        PeptideIdentificationList peptide_ids
        OpenSearchModificationAnalysis_OpenSearchAnalysisResult modification_analysis
        bool is_open_search

        PeptideSearchEngineFIAlgorithm_SearchResult() except + nogil
        PeptideSearchEngineFIAlgorithm_SearchResult(PeptideSearchEngineFIAlgorithm_SearchResult&) except + nogil