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

        # void writePeptideTable(libcpp_vector[ PeptideEntry ] & peptides, libcpp_vector[ size_t ] & reindexed_peptides, libcpp_vector[ PeptideIdentification ] & identifications, String & output_file)
        # void writePeptideTable(libcpp_vector[ PeptideEntry ] & peptides, libcpp_vector[ size_t ] & reindexed_peptides, ConsensusMap & consensus, String & output_file)
        # void writeProteinTable(libcpp_vector[ ProteinEntry ] & proteins, libcpp_vector[ size_t ] & reindexed_proteins, String & output_file)
        # void writeProteinGroups(libcpp_vector[ ISDGroup ] & isd_groups, libcpp_vector[ MSDGroup ] & msd_groups, String & output_file)
        # void countTargetDecoy(libcpp_vector[ MSDGroup ] & msd_groups, ConsensusMap & consensus)
        # void countTargetDecoy(libcpp_vector[ MSDGroup ] & msd_groups, libcpp_vector[ PeptideIdentification ] & peptide_nodes)
        void clearResult()

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":

    cdef cppclass ResolverResult:
        ResolverResult() nogil except +
        ResolverResult(ResolverResult) nogil except +

