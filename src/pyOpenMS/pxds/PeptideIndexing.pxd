from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from FASTAFile cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideIndexing.h>" namespace "OpenMS":
    
    cdef cppclass PeptideIndexing(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        PeptideIndexing() nogil except +
        PeptideIndexing(PeptideIndexing) nogil except + #wrap-ignore

        PeptideIndexing_ExitCodes run(libcpp_vector[ FASTAEntry ] & proteins,
                                      libcpp_vector[ ProteinIdentification ] & prot_ids,
                                      libcpp_vector[ PeptideIdentification ] & pep_ids) nogil except +

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideIndexing.h>" namespace "OpenMS::PeptideIndexing":
    cdef enum PeptideIndexing_ExitCodes "OpenMS::PeptideIndexing::ExitCodes":
        #wrap-attach:
        #    PeptideIndexing
        EXECUTION_OK
        DATABASE_EMPTY
        PEPTIDE_IDS_EMPTY
        ILLEGAL_PARAMETERS
        UNEXPECTED_RESULT

