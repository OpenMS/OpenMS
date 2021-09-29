from Types cimport *
from libcpp.vector cimport vector as libcpp_vector
from DefaultParamHandler cimport *
from ProgressLogger cimport *
from FASTAFile cimport *

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideIndexing.h>" namespace "OpenMS":
    
    cdef cppclass PeptideIndexing(DefaultParamHandler) :
        # wrap-inherits:
        #  DefaultParamHandler
        # wrap-doc:
        #   Refreshes the protein references for all peptide hits in a vector of PeptideIdentifications and adds target/decoy information
        #   -----
        #   All peptide and protein hits are annotated with target/decoy information, using the meta value "target_decoy". For proteins the possible values are "target" and "decoy", 
        #   depending on whether the protein accession contains the decoy pattern (parameter `decoy_string`) as a suffix or prefix, respectively (see parameter `prefix`).
        #   For peptides, the possible values are "target", "decoy" and "target+decoy", depending on whether the peptide sequence is found only in target proteins,
        #   only in decoy proteins, or in both. The target/decoy information is crucial for the @ref TOPP_FalseDiscoveryRate tool.
        #   (For FDR calculations, "target+decoy" peptide hits count as target hits.)

        PeptideIndexing() nogil except +
        PeptideIndexing(PeptideIndexing &) nogil except + # compiler

        PeptideIndexing_ExitCodes run(libcpp_vector[ FASTAEntry ] & proteins,
                                      libcpp_vector[ ProteinIdentification ] & prot_ids,
                                      libcpp_vector[ PeptideIdentification ] & pep_ids) nogil except +

        String getDecoyString() nogil except +
        bool isPrefix() nogil except +

cdef extern from "<OpenMS/ANALYSIS/ID/PeptideIndexing.h>" namespace "OpenMS::PeptideIndexing":
    cdef enum PeptideIndexing_ExitCodes "OpenMS::PeptideIndexing::ExitCodes":
        #wrap-attach:
        #    PeptideIndexing
        EXECUTION_OK
        DATABASE_EMPTY
        PEPTIDE_IDS_EMPTY
        ILLEGAL_PARAMETERS
        UNEXPECTED_RESULT
        DECOYSTRING_EMPTY

