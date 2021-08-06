from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from FASTAFile cimport *
from PreprocessedPairSpectra cimport *
from CrossLinkSpectrumMatch cimport *
from MSExperiment cimport *
from ConsensusMap cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OpenPepXLAlgorithm.h>" namespace "OpenMS":

    cdef cppclass OpenPepXLAlgorithm(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        OpenPepXLAlgorithm() nogil except +
        OpenPepXLAlgorithm(OpenPepXLAlgorithm &) nogil except + # compiler

        OpenPepXLAlgorithm_ExitCodes run(MSExperiment& unprocessed_spectra,
                                         ConsensusMap& cfeatures,
                                         libcpp_vector[ FASTAEntry ]& fasta_db,
                                         libcpp_vector[ ProteinIdentification ]& protein_ids,
                                         libcpp_vector[ PeptideIdentification ]& peptide_ids,
                                         OPXL_PreprocessedPairSpectra& preprocessed_pair_spectra,
                                         libcpp_vector[libcpp_pair[ size_t, size_t ] ]& spectrum_pairs,
                                         libcpp_vector[ libcpp_vector[ CrossLinkSpectrumMatch ] ]& all_top_csms,
                                         MSExperiment& spectra) nogil except +
            # wrap-doc:
                #   Performs the main function of this class, the search for cross-linked peptides
                #   -----
                #   :param unprocessed_spectra: The input PeakMap of experimental spectra
                #   :param cfeatures: The input cfeatures
                #   :param fasta_db: The protein database containing targets and decoys
                #   :param protein_ids: A result vector containing search settings. Should contain one PeptideIdentification
                #   :param peptide_ids: A result vector containing cross-link spectrum matches as PeptideIdentifications and PeptideHits. Should be empty
                #   :param preprocessed_pair_spectra: A result structure containing linear and cross-linked ion spectra. Will be overwritten. This is only necessary for writing out xQuest type spectrum files
                #   :param spectrum_pairs: A result vector containing paired spectra indices. Should be empty. This is only necessary for writing out xQuest type spectrum files
                #   :param all_top_csms: A result vector containing cross-link spectrum matches as CrossLinkSpectrumMatches. Should be empty. This is only necessary for writing out xQuest type spectrum files
                #   :param spectra: A result vector containing the input spectra after preprocessing and filtering. Should be empty. This is only necessary for writing out xQuest type spectrum files

cdef extern from "<OpenMS/ANALYSIS/XLMS/OpenPepXLAlgorithm.h>" namespace "OpenMS::OpenPepXLAlgorithm":
    cdef enum OpenPepXLAlgorithm_ExitCodes "OpenMS::OpenPepXLAlgorithm::ExitCodes":
        #wrap-attach:
        #    OpenPepXLAlgorithm
        EXECUTION_OK
        ILLEGAL_PARAMETERS
        UNEXPECTED_RESULT
        INCOMPATIBLE_INPUT_DATA
