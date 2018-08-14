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
        OpenPepXLAlgorithm(OpenPepXLAlgorithm) nogil except + #wrap-ignore

        OpenPepXLAlgorithm_ExitCodes run(MSExperiment& unprocessed_spectra,
                                         ConsensusMap& cfeatures,
                                         libcpp_vector[ FASTAEntry ]& fasta_db,
                                         libcpp_vector[ ProteinIdentification ]& protein_ids,
                                         libcpp_vector[ PeptideIdentification ]& peptide_ids,
                                         OPXL_PreprocessedPairSpectra& preprocessed_pair_spectra,
                                         libcpp_vector[libcpp_pair[ size_t, size_t ] ]& spectrum_pairs,
                                         libcpp_vector[ libcpp_vector[ CrossLinkSpectrumMatch ] ]& all_top_csms,
                                         MSExperiment& spectra) nogil except +

cdef extern from "<OpenMS/ANALYSIS/XLMS/OpenPepXLAlgorithm.h>" namespace "OpenMS::OpenPepXLAlgorithm":
    cdef enum OpenPepXLAlgorithm_ExitCodes "OpenMS::OpenPepXLAlgorithm::ExitCodes":
        #wrap-attach:
        #    OpenPepXLAlgorithm
        EXECUTION_OK
        ILLEGAL_PARAMETERS
        UNEXPECTED_RESULT
