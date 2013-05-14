from libcpp.vector cimport vector as libcpp_vector

from DefaultParamHandler cimport *

from ProteinIdentification cimport *
from ConsensusMap cimport *
from FASTAFile cimport *

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS":
    cdef cppclass ProteinResolver(DefaultParamHandler):
        # wrap-inherits:
        #    DefaultParamHandler

        ProteinResolver() nogil except +
        ProteinResolver(ProteinResolver) nogil except +

        void resolveConsensus(ConsensusMap & consensus) nogil except +
        void resolveID(libcpp_vector[PeptideIdentification] & peptide_identifications) nogil except +
        void setProteinData(libcpp_vector[FASTAEntry] & protein_data) nogil except +
        libcpp_vector[ResolverResult] getResults() nogil except +

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":

    cdef cppclass ResolverResult:
        ResolverResult() nogil except +
        ResolverResult(ResolverResult) nogil except +

