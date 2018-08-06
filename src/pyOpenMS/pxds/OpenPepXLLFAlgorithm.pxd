from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from libcpp.pair cimport pair as libcpp_pair
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from FASTAFile cimport *
from PreprocessedPairSpectra cimport *
from CrossLinkSpectrumMatch cimport *
from MSExperiment cimport *

cdef extern from "<OpenMS/ANALYSIS/XLMS/OpenPepXLLFAlgorithm.h>" namespace "OpenMS":

    cdef cppclass OpenPepXLLFAlgorithm(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        OpenPepXLLFAlgorithm() nogil except +
        OpenPepXLLFAlgorithm(OpenPepXLLFAlgorithm) nogil except + #wrap-ignore

        OpenPepXLLFAlgorithm_ExitCodes run(MSExperiment& unprocessed_spectra,
                                           libcpp_vector[ FASTAEntry ]& fasta_db,
                                           libcpp_vector[ ProteinIdentification ]& protein_ids,
                                           libcpp_vector[ PeptideIdentification ]& peptide_ids,
                                           libcpp_vector[ libcpp_vector[ CrossLinkSpectrumMatch ] ]& all_top_csms,
                                           MSExperiment& spectra) nogil except +

cdef extern from "<OpenMS/ANALYSIS/XLMS/OpenPepXLLFAlgorithm.h>" namespace "OpenMS::OpenPepXLLFAlgorithm":
    cdef enum OpenPepXLLFAlgorithm_ExitCodes "OpenMS::OpenPepXLLFAlgorithm::ExitCodes":
        #wrap-attach:
        #    OpenPepXLLFAlgorithm
        EXECUTION_OK
        ILLEGAL_PARAMETERS
        UNEXPECTED_RESULT
