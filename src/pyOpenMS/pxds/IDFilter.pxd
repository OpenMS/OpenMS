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

cdef extern from "<OpenMS/PROCESSING/ID/IDFilter.h>" namespace "OpenMS":

    cdef cppclass IDFilter:
         # wrap-doc:
                #  Finds the best-scoring hit in a vector of peptide or protein identifications\n
                #  
                #  This class provides functions for filtering collections of peptide or protein identifications according to various criteria.
                #  It also contains helper functions and classes (functors that implement predicates) that are used in this context.\n
                #  
                #  The filter functions modify their inputs, rather than creating filtered copies.\n
                #  
                #  Most filters work on the hit level, i.e. they remove peptide or protein hits from peptide or protein identifications (IDs).
                #  A few filters work on the ID level instead, i.e. they remove peptide or protein IDs from vectors thereof.
                #  Independent of this, the inputs for all filter functions are vectors of IDs, because the data most often comes in this form.
                #  This design also allows many helper objects to be set up only once per vector, rather than once per ID.\n
                #  
                #  The filter functions for vectors of peptide/protein IDs do not include clean-up steps (e.g. removal of IDs without hits, reassignment of hit ranks, ...).
                #  They only carry out their specific filtering operations.
                #  This is so filters can be chained without having to repeat clean-up operations.
                #  The group of clean-up functions provides helpers that are useful to ensure data integrity after filters have been applied, but it is up to the individual developer to use them when necessary.\n
                #  
                #  The filter functions for MS/MS experiments do include clean-up steps, because they filter peptide and protein IDs in conjunction and potential contradictions between the two must be eliminated.

        IDFilter() except + nogil  
        IDFilter(IDFilter &) except + nogil  #compiler

        Size countHits(libcpp_vector[PeptideIdentification] identifications) except + nogil  # wrap-doc:Returns the total number of peptide hits in a vector of peptide identifications
        Size countHits(libcpp_vector[ProteinIdentification] identifications) except + nogil  # wrap-doc:Returns the total number of protein hits in a vector of protein identifications

        bool getBestHit(libcpp_vector[PeptideIdentification] identifications, bool assume_sorted, PeptideHit& best_hit) except + nogil 
            # wrap-doc:
                #  Finds the best-scoring hit in a vector of peptide or protein identifications\n
                #  
                #  If there are several hits with the best score, the first one is taken
                #  
                #  
                #  :param identifications: Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
                #  :param assume_sorted: Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
                #  :param best_hit: Contains the best hit if successful in a vector of peptide identifications
                #  :return: true if a hit was present, false otherwise

        bool getBestHit(libcpp_vector[ProteinIdentification] identifications, bool assume_sorted, ProteinHit& best_hit) except + nogil 
            # wrap-doc:
                #  Finds the best-scoring hit in a vector of peptide or protein identifications
                #  
                #  If there are several hits with the best score, the first one is taken
                #  
                #  
                #  :param identifications: Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
                #  :param assume_sorted: Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
                #  :param best_hit: Contains the best hit if successful in a vector of protein identifications
                #  :return: true if a hit was present, false otherwise

        void extractPeptideSequences(libcpp_vector[PeptideIdentification]& peptides, libcpp_set[String]& sequences, bool ignore_mods) except + nogil 
            # wrap-doc:
                #  Extracts all unique peptide sequences from a list of peptide IDs
                #  
                #  
                #  :param peptides:
                #  :param ignore_mods: Boolean operator default to false in case of any modifications in sequences during extraction
                #  :return: Sequences

        void updateHitRanks(libcpp_vector[PeptideIdentification]& identifications) except + nogil  # wrap-doc:Updates the hit ranks on all peptide or protein IDs
        void updateHitRanks(libcpp_vector[ProteinIdentification]& identifications) except + nogil  # wrap-doc:Updates the hit ranks on all peptide or protein IDs

        void removeUnreferencedProteins(libcpp_vector[ProteinIdentification]& proteins, libcpp_vector[PeptideIdentification]& peptides) except + nogil  # wrap-doc:Removes protein hits from the protein IDs in a 'cmap' that are not referenced by a peptide in the features or if requested in the unassigned peptide list

        void updateProteinReferences(libcpp_vector[PeptideIdentification]& peptides, libcpp_vector[ProteinIdentification]& proteins, bool remove_peptides_without_reference) except + nogil  # wrap-doc:Removes references to missing proteins. Only PeptideEvidence entries that reference protein hits in 'proteins' are kept in the peptide hits

        bool updateProteinGroups(libcpp_vector[ProteinGroup]& groups, libcpp_vector[ProteinHit]& hits) except + nogil 
            # wrap-doc:
                #  Update protein groups after protein hits were filtered
                #  
                #  
                #  :param groups: Input/output protein groups
                #  :param hits: Available protein hits (all others are removed from the groups)
                #  :return: Returns whether the groups are still valid (which is the case if only whole groups, if any, were removed)

        void removeEmptyIdentifications(libcpp_vector[PeptideIdentification]& ids) except + nogil  # wrap-doc:Removes peptide or protein identifications that have no hits in them
        void removeEmptyIdentifications(libcpp_vector[ProteinIdentification]& ids) except + nogil  # wrap-doc:Removes peptide or protein identifications that have no hits in them

        void filterHitsByScore(libcpp_vector[PeptideIdentification]& ids, double threshold_score) except + nogil  # wrap-doc:Filters peptide or protein identifications according to the score of the hits. The score orientation has to be set to higherscorebetter in each PeptideIdentification. Only peptide/protein hits with a score at least as good as 'threshold_score' are kept
        void filterHitsByScore(libcpp_vector[ProteinIdentification]& ids, double threshold_score) except + nogil  # wrap-doc:Filters peptide or protein identifications according to the score of the hits. The score orientation has to be set to higherscorebetter in each PeptideIdentification/ProteinIdentifiation. Only peptide/protein hits with a score at least as good as 'threshold_score' are kept

        void keepNBestSpectra(libcpp_vector[PeptideIdentification]& peptides, Size n) except + nogil  # wrap-doc:Filter identifications by "N best" PeptideIdentification objects (better PeptideIdentification means better [best] PeptideHit than other)

        void keepNBestHits(libcpp_vector[PeptideIdentification]& ids, Size n) except + nogil  # TODO
        void keepNBestHits(libcpp_vector[ProteinIdentification]& ids, Size n) except + nogil  # TODO

        void filterHitsByRank(libcpp_vector[PeptideIdentification]& ids, Size min_rank, Size max_rank) except + nogil 
            # wrap-doc:
                #  Filters peptide or protein identifications according to the ranking of the hits\n
                #  
                #  The hits between 'min_rank' and 'max_rank' (both inclusive) in each ID are kept
                #  Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1
                #  The ranks are (re-)computed before filtering
                #  'max_rank' is ignored if it is smaller than 'min_rank'
                #  
                #  
                #  Note: There may be several hits with the same rank in a peptide or protein ID (if the scores are the same). This method is useful if a range of higher hits is needed for decoy fairness analysis

        void filterHitsByRank(libcpp_vector[ProteinIdentification]& ids, Size min_rank, Size max_rank) except + nogil 
            # wrap-doc:
                #  Filters peptide or protein identifications according to the ranking of the hits\n
                #  
                #  The hits between 'min_rank' and 'max_rank' (both inclusive) in each ID are kept
                #  Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1
                #  The ranks are (re-)computed before filtering
                #  'max_rank' is ignored if it is smaller than 'min_rank'
                #   
                #  
                #  Note: There may be several hits with the same rank in a peptide or protein ID (if the scores are the same). This method is useful if a range of higher hits is needed for decoy fairness analysis

        void removeDecoyHits(libcpp_vector[PeptideIdentification]& ids) except + nogil  # wrap-doc:Removes hits annotated as decoys from peptide or protein identifications. Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively
        void removeDecoyHits(libcpp_vector[ProteinIdentification]& ids) except + nogil  # wrap-doc:Removes hits annotated as decoys from peptide or protein identifications. Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively

        void removeHitsMatchingProteins(libcpp_vector[PeptideIdentification]& ids, libcpp_set[String] accessions) except + nogil  # wrap-doc:Filters peptide or protein identifications according to the given proteins (negative)
        void removeHitsMatchingProteins(libcpp_vector[ProteinIdentification]& ids, libcpp_set[String] accessions) except + nogil  # wrap-doc:Filters peptide or protein identifications according to the given proteins (negative)

        void keepHitsMatchingProteins(libcpp_vector[PeptideIdentification]& ids, libcpp_set[String] accessions) except + nogil  # wrap-doc:Filters peptide or protein identifications according to the given proteins (positive)
        void keepHitsMatchingProteins(libcpp_vector[ProteinIdentification]& ids, libcpp_set[String] accessions) except + nogil  # wrap-doc:Filters peptide or protein identifications according to the given proteins (positive)

        void keepBestPeptideHits(libcpp_vector[PeptideIdentification]& peptides, bool strict) except + nogil  
            # wrap-doc:
                #  Filters peptide identifications keeping only the single best-scoring hit per ID
                #  
                #  
                #  :param peptides: Input/output
                #  :param strict: If set, keep the best hit only if its score is unique - i.e. ties are not allowed. (Otherwise all hits with the best score is kept.)

        void filterPeptidesByLength(libcpp_vector[PeptideIdentification]& peptides, Size min_length, Size max_length) except + nogil  # wrap-doc:Filters peptide identifications according to peptide sequence length

        void filterPeptidesByCharge(libcpp_vector[PeptideIdentification]& peptides, Size min_charge, Size max_charge) except + nogil  # wrap-doc:Filters peptide identifications according to charge state

        void filterPeptidesByRT(libcpp_vector[PeptideIdentification]& peptides, Size min_rt, Size max_rt) except + nogil  # wrap-doc:Filters peptide identifications by precursor RT, keeping only IDs in the given range

        void filterPeptidesByMZ(libcpp_vector[PeptideIdentification]& peptides, Size min_mz, Size max_mz) except + nogil  # wrap-doc:Filters peptide identifications by precursor m/z, keeping only IDs in the given range

        void filterPeptidesByMZError(libcpp_vector[PeptideIdentification]& peptides, double mass_error, bool unit_ppm) except + nogil  # wrap-doc:Filter peptide identifications according to mass deviation

        void filterPeptidesByRTPredictPValue(libcpp_vector[PeptideIdentification]& peptides, const String& metavalue_key, double threshold) except + nogil 
            # wrap-doc:
                #  Filters peptide identifications according to p-values from RTPredict\n
                #  
                #  Filters the peptide hits by the probability (p-value) of a correct peptide identification having a deviation between observed and predicted RT equal to or greater than allowed
                #  
                #  
                #  :param peptides: Input/output
                #  :param metavalue_key: Name of the meta value that holds the p-value: "predicted_RT_p_value" or "predicted_RT_p_value_first_dim"
                #  :param threshold: P-value threshold

        void removePeptidesWithMatchingModifications(libcpp_vector[PeptideIdentification]& peptides, libcpp_set[String]& modifications) except + nogil  # wrap-doc:Removes all peptide hits that have at least one of the given modifications

        void keepPeptidesWithMatchingModifications(libcpp_vector[PeptideIdentification]& peptides, libcpp_set[String]& modifications) except + nogil  # wrap-doc:Keeps only peptide hits that have at least one of the given modifications

        void removePeptidesWithMatchingSequences(libcpp_vector[PeptideIdentification]& peptides, libcpp_vector[PeptideIdentification]& bad_peptides, bool ignore_mods) except + nogil  # wrap-doc:Removes all peptide hits with a sequence that matches one in 'bad_peptides'

        void keepPeptidesWithMatchingSequences(libcpp_vector[PeptideIdentification]& peptides, libcpp_vector[PeptideIdentification]& bad_peptides, bool ignore_mods) except + nogil  # wrap-doc:Removes all peptide hits with a sequence that does not match one in 'good_peptides'

        void keepUniquePeptidesPerProtein(libcpp_vector[PeptideIdentification]& peptides) except + nogil  # wrap-doc:Removes all peptides that are not annotated as unique for a protein (by PeptideIndexer)

        void removeDuplicatePeptideHits(libcpp_vector[PeptideIdentification]& peptides) except + nogil  # wrap-doc:Removes duplicate peptide hits from each peptide identification, keeping only unique hits (per ID)

        void filterHitsByScore(MSExperiment& experiment, double peptide_threshold_score, double protein_threshold_score) except + nogil  # wrap-doc:Filters an MS/MS experiment according to score thresholds

        void keepNBestHits(MSExperiment& experiment, Size n) except + nogil  # wrap-doc:Filters an MS/MS experiment by keeping the N best peptide hits for every spectrum
        
        void keepBestPerPeptide(libcpp_vector[PeptideIdentification]& peptides, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum) except + nogil  # wrap-doc:Filters PeptideHits from PeptideIdentification by keeping only the best peptide hits for every peptide sequence

        void keepBestPerPeptidePerRun(libcpp_vector[ProteinIdentification]& prot_ids, libcpp_vector[PeptideIdentification]& peptides, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum) except + nogil  # wrap-doc:Filters PeptideHits from PeptideIdentification by keeping only the best peptide hits for every peptide sequence on a per run basis
        
        void keepHitsMatchingProteins(MSExperiment& experiment, libcpp_vector[FASTAEntry]& proteins) except + nogil 

cdef extern from "<OpenMS/PROCESSING/ID/IDFilter.h>" namespace "OpenMS::IDFilter":
    
    cdef cppclass DigestionFilter "OpenMS::IDFilter::DigestionFilter":
        # wrap-attach:
        #   IDFilter
        DigestionFilter() except + nogil  # wrap-ignore
        DigestionFilter(DigestionFilter) except + nogil  #wrap-ignore

        # GetMatchingItems[ PeptideEvidence, FASTAEntry ] accession_resolver_
        ProteaseDigestion digestion_
        bool ignore_missed_cleavages_
        bool methionine_cleavage_

        DigestionFilter(libcpp_vector[ FASTAEntry ] & entries, 
                        ProteaseDigestion & digestion,
                        bool ignore_missed_cleavages,
                        bool methionine_cleavage) except + nogil 

        # bool operator()(PeptideEvidence & evidence) except + nogil 
        void filterPeptideEvidences(libcpp_vector[ PeptideIdentification ] & peptides) except + nogil 
