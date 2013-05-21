from libcpp.vector cimport vector as libcpp_vector
from libcpp.list cimport list as libcpp_list
from DefaultParamHandler cimport *
from ProteinIdentification cimport *
from ConsensusMap cimport *
from FASTAFile cimport *
from String cimport *

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

        # TODO not defined in OpenMS
        # void writePeptideTable(libcpp_vector[ PeptideEntry ] & peptides, libcpp_vector[ size_t ] & reindexed_peptides, libcpp_vector[ PeptideIdentification ] & identifications, String & output_file)
        # void writePeptideTable(libcpp_vector[ PeptideEntry ] & peptides, libcpp_vector[ size_t ] & reindexed_peptides, ConsensusMap & consensus, String & output_file)
        # void writeProteinTable(libcpp_vector[ ProteinEntry ] proteins, libcpp_vector[ size_t ] & reindexed_proteins, String & output_file)
        # void writeProteinGroups(libcpp_vector[ ISDGroup ] & isd_groups, libcpp_vector[ MSDGroup ] & msd_groups, String & output_file)
        void countTargetDecoy(libcpp_vector[ MSDGroup ] & msd_groups, ConsensusMap & consensus)
        void countTargetDecoy(libcpp_vector[ MSDGroup ] & msd_groups, libcpp_vector[ PeptideIdentification ] & peptide_nodes)
        void clearResult()

        # POINTER # PeptideIdentification  getPeptideIdentification(ConsensusMap & consensus, PeptideEntry * peptide)
        # POINTER # PeptideHit  getPeptideHit(ConsensusMap & consensus, PeptideEntry * peptide)
        # POINTER # PeptideIdentification  getPeptideIdentification(libcpp_vector[ PeptideIdentification ] & peptide_nodes, PeptideEntry * peptide)
        # POINTER # PeptideHit  getPeptideHit(libcpp_vector[ PeptideIdentification ] & peptide_nodes, PeptideEntry * peptide)

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass ISDGroup "OpenMS::ProteinResolver::ISDGroup":
        ISDGroup() nogil except +
        ISDGroup(ISDGroup) nogil except + #wrap-ignore
        # NAMESPACE # # POINTER # std::list[ ProteinEntry * ] proteins
        # NAMESPACE # # POINTER # std::list[ PeptideEntry * ] peptides
        Size index
        # libcpp_list[ size_t ] msd_groups # TODO no converter for List

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass MSDGroup "OpenMS::ProteinResolver::MSDGroup":
        MSDGroup() nogil except +
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

    cdef cppclass ResolverResult "OpenMS::ProteinResolver::ResolverResult":
        ResolverResult() nogil except +
        ResolverResult(ResolverResult) nogil except +

        String identifier
        # std::vector<ISDGroup> * isds;
        # std::vector<MSDGroup> * msds;
        # std::vector<ProteinEntry> * protein_entries;
        # std::vector<PeptideEntry> * peptide_entries;
        # std::vector<Size> * reindexed_peptides;
        # std::vector<Size> * reindexed_proteins;
        # enum type  {PeptideIdent, Consensus} input_type;
        # std::vector<PeptideIdentification> * peptide_identification;
        # ConsensusMap * consensus_map;

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass PeptideEntry "OpenMS::ProteinResolver::PeptideEntry":
        PeptideEntry() nogil except +
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
        ProteinEntry() nogil except +
        ProteinEntry(ProteinEntry) nogil except + #wrap-ignore
        # NAMESPACE # # POINTER # std::list[ PeptideEntry * ] peptides
        bool traversed
        # NAMESPACE # # POINTER # FASTAFile::FASTAEntry * fasta_entry
        ProteinEntry_type protein_type
        DoubleReal weight
        Real coverage
        # NAMESPACE # # POINTER # std::list[ ProteinEntry * ] indis
        Size index
        Size msd_group
        Size isd_group
        Size number_of_experimental_peptides

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver::ProteinEntry":
    cdef enum ProteinEntry_type "OpenMS::ProteinResolver::ProteinEntry::type":
        #wrap-attach:
        #    ProteinEntry
        primary
        secondary
        primary_indistinguishable
        secondary_indistinguishable

