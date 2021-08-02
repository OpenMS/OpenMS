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
        ProteinResolver(ProteinResolver &) nogil except +

        void resolveConsensus(ConsensusMap & consensus) nogil except +
        void resolveID(libcpp_vector[PeptideIdentification] & peptide_identifications) nogil except +
        void setProteinData(libcpp_vector[FASTAEntry] & protein_data) nogil except +
        libcpp_vector[ResolverResult] getResults() nogil except +

        void countTargetDecoy(libcpp_vector[ MSDGroup ] & msd_groups, ConsensusMap & consensus) nogil except +
        void countTargetDecoy(libcpp_vector[ MSDGroup ] & msd_groups, libcpp_vector[ PeptideIdentification ] & peptide_nodes) nogil except +
        void clearResult() nogil except +

        PeptideIdentification getPeptideIdentification(ConsensusMap & consensus, PeptideEntry * peptide) nogil except +
        PeptideHit getPeptideHit(ConsensusMap & consensus, PeptideEntry * peptide) nogil except +
        PeptideIdentification getPeptideIdentification(libcpp_vector[ PeptideIdentification ] & peptide_nodes, PeptideEntry * peptide) nogil except +
        PeptideHit getPeptideHit(libcpp_vector[ PeptideIdentification ] & peptide_nodes, PeptideEntry * peptide) nogil except +

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass ISDGroup "OpenMS::ProteinResolver::ISDGroup":
        ISDGroup() nogil except + # compiler
        ISDGroup(ISDGroup &) nogil except + # compiler

        # NAMESPACE # # POINTER # std::list[ ProteinEntry * ] proteins
        # NAMESPACE # # POINTER # std::list[ PeptideEntry * ] peptides
        Size index
        # libcpp_list[ size_t ] msd_groups # TODO no converter for List

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass MSDGroup "OpenMS::ProteinResolver::MSDGroup":
        MSDGroup() nogil except + # compiler
        MSDGroup(MSDGroup &) nogil except + # compiler

        # NAMESPACE # # POINTER # std::list[ ProteinEntry * ] proteins
        # NAMESPACE # # POINTER # std::list[ PeptideEntry * ] peptides
        Size index
        ISDGroup * isd_group
        Size number_of_decoy
        Size number_of_target
        Size number_of_target_plus_decoy
        float intensity

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":

    cdef cppclass ResolverResult "OpenMS::ProteinResolver::ResolverResult":
        ResolverResult() nogil except + # compiler
        ResolverResult(ResolverResult &) nogil except + # compiler

        String identifier
        libcpp_vector[ISDGroup] * isds
        libcpp_vector[MSDGroup] * msds
        libcpp_vector[ProteinEntry] * protein_entries
        libcpp_vector[PeptideEntry] * peptide_entries
        # libcpp_vector[size_t] * reindexed_peptides
        # libcpp_vector[size_t] * reindexed_proteins
        ProteinResolverResult_Type input_type
        libcpp_vector[PeptideIdentification] * peptide_identification
        ConsensusMap * consensus_map

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass PeptideEntry "OpenMS::ProteinResolver::PeptideEntry":
        PeptideEntry() nogil except + # compiler
        PeptideEntry(PeptideEntry &) nogil except + # compiler

        # NAMESPACE # # POINTER # std::list[ ProteinEntry * ] proteins
        bool traversed
        String sequence
        Size peptide_identification
        Size peptide_hit
        Size index
        Size msd_group
        Size isd_group
        bool experimental
        float intensity
        String origin

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver::ResolverResult":
    
    cdef enum ProteinResolverResult_Type "OpenMS::ProteinResolver::ResolverResult::type":
        PeptideIdent
        Consensus

cdef extern from "<OpenMS/ANALYSIS/QUANTITATION/ProteinResolver.h>" namespace "OpenMS::ProteinResolver":
    
    cdef cppclass ProteinEntry "OpenMS::ProteinResolver::ProteinEntry":
        ProteinEntry() nogil except + # compiler
        ProteinEntry(ProteinEntry &) nogil except + # compiler

        # NAMESPACE # # POINTER # std::list[ PeptideEntry * ] peptides
        bool traversed
        FASTAEntry * fasta_entry
        ProteinEntry_type protein_type
        double weight
        float coverage
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

