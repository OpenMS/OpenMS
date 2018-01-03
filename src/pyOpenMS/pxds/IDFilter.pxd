from libcpp.vector cimport vector as libcpp_vector
from libcpp.set cimport set as libcpp_set
from libcpp cimport bool

from String cimport *
from DefaultParamHandler cimport *

from PeptideIdentification cimport *
from ProteinIdentification cimport *
from FASTAFile cimport *
from ProteaseDigestion cimport *

from MSExperiment cimport *
from MSSpectrum cimport *
from Peak1D cimport *
from ChromatogramPeak cimport *

cdef extern from "<OpenMS/FILTERING/ID/IDFilter.h>" namespace "OpenMS":

    cdef cppclass IDFilter:

        IDFilter()           nogil except +
        IDFilter(IDFilter)   nogil except + # wrap-ignore

        Size countHits(libcpp_vector[PeptideIdentification] identifications) nogil except +
        Size countHits(libcpp_vector[ProteinIdentification] identifications) nogil except +

        bool getBestHit(libcpp_vector[PeptideIdentification] identifications, bool assume_sorted, PeptideHit& best_hit) nogil except +
        bool getBestHit(libcpp_vector[ProteinIdentification] identifications, bool assume_sorted, ProteinHit& best_hit) nogil except +

        void extractPeptideSequences(libcpp_vector[PeptideIdentification]& peptides, libcpp_set[String]& sequences, bool ignore_mods) nogil except +

        void updateHitRanks(libcpp_vector[PeptideIdentification]& identifications) nogil except +
        void updateHitRanks(libcpp_vector[ProteinIdentification]& identifications) nogil except +

        void removeUnreferencedProteins(libcpp_vector[ProteinIdentification]& proteins, libcpp_vector[PeptideIdentification]& peptides) nogil except +

        void updateProteinReferences(libcpp_vector[PeptideIdentification]& peptides, libcpp_vector[ProteinIdentification]& proteins, bool remove_peptides_without_reference) nogil except +

        bool updateProteinGroups(libcpp_vector[ProteinGroup]& groups, libcpp_vector[ProteinHit]& hits) nogil except +

        void removeEmptyIdentifications(libcpp_vector[PeptideIdentification]& ids) nogil except +
        void removeEmptyIdentifications(libcpp_vector[ProteinIdentification]& ids) nogil except +

        void filterHitsByScore(libcpp_vector[PeptideIdentification]& ids, double threshold_score) nogil except +
        void filterHitsByScore(libcpp_vector[ProteinIdentification]& ids, double threshold_score) nogil except +

        void filterHitsBySignificance(libcpp_vector[PeptideIdentification]& ids, double threshold_fraction) nogil except +
        void filterHitsBySignificance(libcpp_vector[ProteinIdentification]& ids, double threshold_fraction) nogil except +

        void keepNBestHits(libcpp_vector[PeptideIdentification]& ids, Size n) nogil except +
        void keepNBestHits(libcpp_vector[ProteinIdentification]& ids, Size n) nogil except +

        void filterHitsByRank(libcpp_vector[PeptideIdentification]& ids, Size min_rank, Size max_rank) nogil except +
        void filterHitsByRank(libcpp_vector[ProteinIdentification]& ids, Size min_rank, Size max_rank) nogil except +

        void removeDecoyHits(libcpp_vector[PeptideIdentification]& ids) nogil except +
        void removeDecoyHits(libcpp_vector[ProteinIdentification]& ids) nogil except +

        void removeHitsMatchingProteins(libcpp_vector[PeptideIdentification]& ids, libcpp_set[String] accessions) nogil except +
        void removeHitsMatchingProteins(libcpp_vector[ProteinIdentification]& ids, libcpp_set[String] accessions) nogil except +

        void keepHitsMatchingProteins(libcpp_vector[PeptideIdentification]& ids, libcpp_set[String] accessions) nogil except +
        void keepHitsMatchingProteins(libcpp_vector[ProteinIdentification]& ids, libcpp_set[String] accessions) nogil except +

        void keepBestPeptideHits(libcpp_vector[PeptideIdentification]& peptides, bool strict) nogil except +

        void filterPeptidesByLength(libcpp_vector[PeptideIdentification]& peptides, Size min_length, Size max_length) nogil except +

        void filterPeptidesByCharge(libcpp_vector[PeptideIdentification]& peptides, Size min_charge, Size max_charge) nogil except +

        void filterPeptidesByRT(libcpp_vector[PeptideIdentification]& peptides, Size min_rt, Size max_rt) nogil except +

        void filterPeptidesByMZ(libcpp_vector[PeptideIdentification]& peptides, Size min_mz, Size max_mz) nogil except +

        void filterPeptidesByMZError(libcpp_vector[PeptideIdentification]& peptides, double mass_error, bool unit_ppm) nogil except +

        void filterPeptidesByRTPredictPValue(libcpp_vector[PeptideIdentification]& peptides, String& metavalue_key, double threshold) nogil except +

        void removePeptidesWithMatchingModifications(libcpp_vector[PeptideIdentification]& peptides, libcpp_set[String]& modifications) nogil except +

        void keepPeptidesWithMatchingModifications(libcpp_vector[PeptideIdentification]& peptides, libcpp_set[String]& modifications) nogil except +

        void removePeptidesWithMatchingSequences(libcpp_vector[PeptideIdentification]& peptides, libcpp_vector[PeptideIdentification]& bad_peptides, bool ignore_mods) nogil except +

        void keepPeptidesWithMatchingSequences(libcpp_vector[PeptideIdentification]& peptides, libcpp_vector[PeptideIdentification]& bad_peptides, bool ignore_mods) nogil except +

        void keepUniquePeptidesPerProtein(libcpp_vector[PeptideIdentification]& peptides) nogil except +

        void removeDuplicatePeptideHits(libcpp_vector[PeptideIdentification]& peptides) nogil except +

        void filterHitsByScore(MSExperiment& experiment, double peptide_threshold_score, double protein_threshold_score) nogil except +

        void filterHitsBySignificance(MSExperiment& experiment, double peptide_threshold_fraction, double protein_threshold_fraction) nogil except +

        void keepNBestHits(MSExperiment& experiment, Size n) nogil except +

        void keepHitsMatchingProteins(MSExperiment& experiment, libcpp_vector[FASTAEntry]& proteins) nogil except +


cdef extern from "<OpenMS/FILTERING/ID/IDFilter.h>" namespace "OpenMS::IDFilter":
    
    cdef cppclass DigestionFilter "OpenMS::IDFilter::DigestionFilter":
        # wrap-attach:
        #    IDFilter
        DigestionFilter() nogil except + # wrap-ignore
        DigestionFilter(DigestionFilter) nogil except + #wrap-ignore

        # GetMatchingItems[ PeptideEvidence, FASTAEntry ] accession_resolver_
        ProteaseDigestion digestion_
        bool ignore_missed_cleavages_
        bool methionine_cleavage_

        DigestionFilter(libcpp_vector[ FASTAEntry ] & entries, 
                        ProteaseDigestion & digestion,
                        bool ignore_missed_cleavages,
                        bool methionine_cleavage) nogil except +

        # bool operator()(PeptideEvidence & evidence) nogil except +
        void filterPeptideEvidences(libcpp_vector[ PeptideIdentification ] & peptides) nogil except +

