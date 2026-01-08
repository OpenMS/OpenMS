from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ProteinIdentification cimport *
from PeptideIdentification cimport *
from String cimport *
from ProgressLogger cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h>" namespace "OpenMS":
 
    cdef cppclass PeptideSearchEngineFIAlgorithm(DefaultParamHandler, ProgressLogger):
        # wrap-inherits:
        #   DefaultParamHandler
        #   ProgressLogger
        # wrap-doc:
        #   Fragment-index-based peptide database search algorithm (experimental).
        #   
        #   Provides a self-contained search engine that matches MS/MS spectra against a protein
        #   database using an FI (Fragment Index). Typical usage:
        #   - Configure parameters via DefaultParamHandler (mass tolerances, enzyme, charges, etc.)
        #   - Call search() with an input mzML file and a FASTA database to populate identification
        #     outputs (ProteinIdentification and PeptideIdentificationList)
        #   - Intended for educational/prototyping use and to demonstrate FI-backed searching

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

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideSearchEngineFIAlgorithm.h>" namespace "OpenMS::PeptideSearchEngineFIAlgorithm":
    cdef enum PeptideSearchEngineFIAlgorithm_ExitCodes "OpenMS::PeptideSearchEngineFIAlgorithm::ExitCodes":
        #wrap-attach:
        #   PeptideSearchEngineFIAlgorithm
        EXECUTION_OK
        INPUT_FILE_EMPTY
        UNEXPECTED_RESULT
        UNKNOWN_ERROR
        ILLEGAL_PARAMETERS