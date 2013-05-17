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

        # TODO
        # void writePeptideTable(libcpp_vector[ PeptideEntry ] & peptides, libcpp_vector[ size_t ] & reindexed_peptides, libcpp_vector[ PeptideIdentification ] & identifications, String & output_file)
        # void writePeptideTable(libcpp_vector[ PeptideEntry ] & peptides, libcpp_vector[ size_t ] & reindexed_peptides, ConsensusMap & consensus, String & output_file)
        # void writeProteinTable(libcpp_vector[ ProteinEntry ] & proteins, libcpp_vector[ size_t ] & reindexed_proteins, String & output_file)
        # void writeProteinGroups(libcpp_vector[ ISDGroup ] & isd_groups, libcpp_vector[ MSDGroup ] & msd_groups, String & output_file)
        # void countTargetDecoy(libcpp_vector[ MSDGroup ] & msd_groups, ConsensusMap & consensus)
        # void countTargetDecoy(libcpp_vector[ MSDGroup ] & msd_groups, libcpp_vector[ PeptideIdentification ] & peptide_nodes)
        void clearResult()

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass ISDGroup "OpenMS::ProteinResolver::ISDGroup":
        ISDGroup(ISDGroup) nogil except + #wrap-ignore
        # NAMESPACE # # POINTER # std::list[ ProteinEntry * ] proteins
        # NAMESPACE # # POINTER # std::list[ PeptideEntry * ] peptides
        Size index
        # NAMESPACE # std::list[ size_t ] msd_groups

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass MSDGroup "OpenMS::ProteinResolver::MSDGroup":
        MSDGroup(MSDGroup) nogil except + #wrap-ignore
        # NAMESPACE # # POINTER # std::list[ ProteinEntry * ] proteins
        # NAMESPACE # # POINTER # std::list[ PeptideEntry * ] peptides
        Size index
        # POINTER # ISDGroup * isd_group
        Size number_of_decoy
        Size number_of_target
        Size number_of_target_plus_decoy
        Real intensity

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":

    cdef cppclass ResolverResult:
        ResolverResult() nogil except +
        ResolverResult(ResolverResult) nogil except +

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass PeptideEntry "OpenMS::ProteinResolver::PeptideEntry":
        PeptideEntry(PeptideEntry) nogil except + #wrap-ignore
        # NAMESPACE # # POINTER # std::list[ ProteinEntry * ] proteins
        bool traversed
        String sequence
        Size peptide_identification
        Size peptide_hit
        Size index
        Size msd_group
        Size isd_group
        bool experimental
        Real intensity
        String origin

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver::ResolverResult":
    
    cdef enum ProteinResolverResult_Type "OpenMS::ProteinResolver::ResolverResult::type":
        PeptideIdent
        Consensus

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass ProteinEntry "OpenMS::ProteinResolver::ProteinEntry":
        ProteinEntry(ProteinEntry) nogil except + #wrap-ignore
        # NAMESPACE # # POINTER # std::list[ PeptideEntry * ] peptides
        bool traversed
        # NAMESPACE # # POINTER # FASTAFile::FASTAEntry * fasta_entry
        # NAMESPACE # enum OpenMS::ProteinResolver::ProteinEntry::type protein_type
        DoubleReal weight
        Real coverage
        # NAMESPACE # # POINTER # std::list[ ProteinEntry * ] indis
        Size index
        Size msd_group
        Size isd_group
        Size number_of_experimental_peptides

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver::ProteinEntry":
    cdef enum type "OpenMS::ProteinResolver::ProteinEntry::type":
        #wrap-attach:
        #    ProteinEntry
        primary
        secondary
        primary_indistinguishable
        secondary_indistinguishable

