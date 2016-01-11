// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Nico Pfeifer, Mathias Walzer, Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_ID_IDFILTER_H
#define OPENMS_FILTERING_ID_IDFILTER_H

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <algorithm>
#include <climits>
#include <functional>
#include <vector>

namespace OpenMS
{
  /**
    @brief Used to filter identifications by different criteria.

    The identifications are filtered by significance thresholds and
    by sequences. The filtering by significance thresholds looks for the
    best ProteinIdentification that fulfills the significance threshold criterion.
    score > significance-threshold * significance_fraction.
    The filtering by sequences looks for the best ProteinIdentification that
    is contained in one of the protein sequences.

    TODO: fix design of filter functions. There will be an error e.g. if input and output points to the same PeptideIdentification.
  */
  class OPENMS_DLLAPI IDFilter
  {
public:

    /// Constructor
    IDFilter();

    /// Destructor
    virtual ~IDFilter();

    /// Predicate for peptide hits: does the sequence have at least the given length?
    struct HasMinPeptideLength;

    /// Predicate for peptide hits: does the charge have at least the given value?
    struct HasMinCharge;
    
    /// Predicate for peptide hits: is the m/z error at least this low?
    struct HasLowMZError;

    /**
       @brief Predicate for peptide hits: given a list of modifications, do any occur in the sequence?

       If the list of modifications is empty, return true if the sequence is modified at all.
    */
    struct HasMatchingModification;

    /// Predicate for peptide hits: given a list of sequences, does one match?
    struct HasMatchingSequence;

    /// Predicate for peptide hits: is the list of peptide evidences empty?
    struct HasNoEvidence;

    /// Predicate for peptide IDs: is the retention time in the given range?
    struct HasRTInRange;

    /// Predicate for peptide IDs: is the retention time in the given range?
    struct HasMZInRange;

    /// Predicate for peptide or protein hits: is the score at least as good as the given value?
    template <class HitType>
    struct HasGoodScore
    {
      typedef HitType argument_type; // for use as a predicate
      
      double score;
      bool higher_score_better;

      HasGoodScore(double score, bool higher_score_better):
        score(score), higher_score_better(higher_score_better)
      {}

      bool operator()(const HitType& hit) const
      {
        if (higher_score_better)
        {
          return hit.getScore() >= score;
        }
        return hit.getScore() <= score;
      }
    };


    /**
       @brief Predicate for peptide or protein hits: is a meta value with given key and value set?

       If the value is empty (DataValue::EMPTY), only the existence of a meta value with the given key is checked.
    */
    template <class HitType>
    struct HasMetaValue
    {
      typedef HitType argument_type; // for use as a predicate

      String key;
      DataValue value;

      HasMetaValue(const String& key, const DataValue& value):
        key(key), value(value)
      {}

      bool operator()(const HitType& hit) const
      {
        DataValue found = hit.getMetaValue(key);
        if (found.isEmpty()) return false; // meta value "key" not set
        if (value.isEmpty()) return true; // "key" is set, value doesn't matter
        return found == value;
      }
    };


    /**
       @brief Predicate for peptide or protein hits: does a meta value have at most the given value?
    */
    template <class HitType>
    struct HasMaxMetaValue
    {
      typedef HitType argument_type; // for use as a predicate

      String key;
      double value;

      HasMaxMetaValue(const String& key, const double& value):
        key(key), value(value)
      {}

      bool operator()(const HitType& hit) const
      {
        DataValue found = hit.getMetaValue(key);
        if (found.isEmpty()) return false; // meta value "key" not set
        return double(found) <= value;
      }
    };


    /// Predicate for peptide or protein hits: is this a decoy hit?
    template <class HitType>
    struct HasDecoyAnnotation
    {
      typedef HitType argument_type; // for use as a predicate

      struct HasMetaValue<HitType> target_decoy, is_decoy;

      HasDecoyAnnotation():
        target_decoy("target_decoy", "decoy"), is_decoy("isDecoy", "true")
      {}

      bool operator()(const HitType& hit) const
      {
        // @TODO: this could be done slightly more efficiently by returning
        // false if the "target_decoy" meta value is "target" or "target+decoy",
        // without checking for an "isDecoy" meta value in that case
        return target_decoy(hit) || is_decoy(hit);
      }
    };


    /// Predicate for peptide/protein hits and peptide evidence: given a list of protein accessions, do any occur in the annotation(s)?
    template <class HitType>
    struct HasMatchingAccession
    {
      typedef HitType argument_type; // for use as a predicate

      const std::set<String>& accessions;
      
      HasMatchingAccession(const std::set<String>& accessions):
        accessions(accessions)
      {}

      bool operator()(const PeptideHit& hit) const
      {
        std::set<String> present_accessions = hit.extractProteinAccessions();
        for (std::set<String>::iterator it = present_accessions.begin();
             it != present_accessions.end(); ++it)
        {
          if (accessions.count(*it) > 0) return true;
        }
        return false;
      }
      
      bool operator()(const ProteinHit& hit) const
      {
        return (accessions.count(hit.getAccession()) > 0);
      }

      bool operator()(const PeptideEvidence& evidence) const
      {
        return (accessions.count(evidence.getProteinAccession()) > 0);
      }
    };


    /// Predicate for peptide/protein identifications: is the list of hits empty?
    template <class IdentificationType>
    struct HasNoHits
    {
      typedef IdentificationType argument_type; // for use as a predicate

      bool operator()(const IdentificationType& id) const
      {
        return id.getHits().empty();
      }
    };


    /// Remove items that satisfy a condition from a container (e.g. vector)
    template <class Container, class Predicate>
    static void removeMatchingItems(Container& items, const Predicate& pred)
    {
      items.erase(std::remove_if(items.begin(), items.end(), pred),
                  items.end());
    }


    /// Keep items that satisfy a condition in a container (e.g. vector), removing all others
    template <class Container, class Predicate>
    static void keepMatchingItems(Container& items, const Predicate& pred)
    {
      items.erase(std::remove_if(items.begin(), items.end(), std::not1(pred)),
                  items.end());
    }


    /**
       @brief Finds the best-scoring hit in a vector of peptide or protein identifications
       
       @param identifications Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
       @param assume_sorted are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
       
       @except Exception::InvalidValue if the IDs have different score types (i.e. scores cannot be compared)
       
       @return true if a hit was present, false otherwise
    */
    template <class IdentificationType>
    static bool getBestHit(const std::vector<IdentificationType> identifications, bool assume_sorted, typename IdentificationType::HitType& best_hit)
    {
      if (identifications.empty()) return false;

      typename std::vector<IdentificationType>::const_iterator best_id_it = identifications.end();
      typename std::vector<typename IdentificationType::HitType>::const_iterator best_hit_it;

      for (typename std::vector<IdentificationType>::const_iterator id_it = identifications.begin(); id_it != identifications.end(); ++id_it)
      {
        if (id_it->getHits().empty()) continue;
        
        if (best_id_it == identifications.end()) // no previous "best" hit
        {
          best_id_it = id_it;
          best_hit_it = id_it->getHits().begin();
        }
        else if (best_id_it->getScoreType() != id_it->getScoreType())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Can't compare scores of different types",
                                        best_id_it->getScoreType() + "/" + id_it->getScoreType());
        }

        bool higher_better = best_id_it->isHigherScoreBetter();
        for (typename std::vector<typename IdentificationType::HitType>::const_iterator hit_it = id_it->getHits().begin();
             hit_it != id_it->getHits().end(); ++hit_it)
        {
          if ((higher_better && (hit_it->getScore() > best_hit_it->getScore())) ||
              (!higher_better && (hit_it->getScore() < best_hit_it->getScore())))
          {
            best_hit_it = hit_it;
          }
          if (assume_sorted) break; // only consider the first hit
        }
      }

      if (best_id_it == identifications.end()) return false; // no hits in any IDs

      best_hit = *best_hit_it;
      return true;
    }


    /// Filters peptide or protein identifications by only allowing peptide/protein hits which reach a score above (or below, depending on score orientation) @p threshold_fraction * significance_threshold
    template <class IdentificationType>
    static void filterHitsBySignificance(std::vector<IdentificationType>& ids,
                                         double threshold_fraction = 1.0)
    {
      for (typename std::vector<IdentificationType>::iterator id_it =
             ids.begin(); id_it != ids.end(); ++id_it)
      {
        double threshold_score = (threshold_fraction *
                                  id_it->getSignificanceThreshold());
        struct HasGoodScore<typename IdentificationType::HitType> score_filter(
          threshold_score, id_it->isHigherScoreBetter());
        keepMatchingItems(id_it->getHits(), score_filter);
      }
    }


    /// Remove peptide or protein identifications that have no hits in them
    template <class IdentificationType>
    static void removeEmptyIdentifications(std::vector<IdentificationType>& ids)
    {
      struct HasNoHits<IdentificationType> empty_filter;
      removeMatchingItems(ids, empty_filter);
    }


    /**
      @brief Filters peptide or protein identifications according to the score of the hits.

      Only peptide/protein hits with a score at least as good as @p threshold_score are kept. Score orientation (are higher scores better?) is taken into account.
    */
    template <class IdentificationType>
    static void filterHitsByScore(std::vector<IdentificationType>& ids,
                                  double threshold_score)
    {
      for (typename std::vector<IdentificationType>::iterator id_it = 
             ids.begin(); id_it != ids.end(); ++id_it)
      {
        struct HasGoodScore<typename IdentificationType::HitType> score_filter(
          threshold_score, id_it->isHigherScoreBetter());
        keepMatchingItems(id_it->getHits(), score_filter);
      }
    }


    /**
      @brief Filters peptide or protein identifications according to the score of the hits, keeping the @p n best hits per ID.

      The score orientation (are higher scores better?) is taken into account.
    */
    template <class IdentificationType>
    static void keepNBestHits(std::vector<IdentificationType>& ids, Size n)
    {
      for (typename std::vector<IdentificationType>::iterator id_it =
             ids.begin(); id_it != ids.end(); ++id_it)
      {
        id_it->sort();
        if (n < id_it->getHits().size()) id_it->getHits().resize(n);
      }
    }


    /**
       @brief Filters peptide or protein identifications according to the ranking of the hits.

       The hits between @p min_rank and @p max_rank (both inclusive) in each ID are kept. Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1. The hits are ranked according to their scores, thus the initial values of the "rank" attribute (PeptideHit::getRank or ProteinHit::getRank) are ignored.

       This method is useful if a range of higher hits are used for decoy fairness analysis.

       @note The ranks of the hits may be invalidated.
    */
    template <class IdentificationType>
    static void filterHitsByRank(std::vector<IdentificationType>& ids,
                                 Size min_rank, Size max_rank)
    {
      for (typename std::vector<IdentificationType>::iterator id_it =
             ids.begin(); id_it != ids.end(); ++id_it)
      {
        id_it->sort();
        std::vector<typename IdentificationType::HitType>& hits = 
          id_it->getHits();
        if ((max_rank > 0) && (max_rank < hits.size())) hits.resize(max_rank);
        if (min_rank > 1) hits.erase(hits.begin(), hits.begin() + min_rank - 1);
      }
    }
    
    
    /**
       @brief Removes hits annotated as decoys from peptide or protein identifications.
     
       Checks for meta values named "target_decoy" and "isDecoy", and removes protein/peptide hits if the values are "decoy" and "true", respectively.

       @note The ranks of the hits may be invalidated.
    */
    template <class IdentificationType>
    static void removeDecoyHits(std::vector<IdentificationType>& ids)
    {
      struct HasDecoyAnnotation<typename IdentificationType::HitType> 
        decoy_filter;
      for (typename std::vector<IdentificationType>::iterator id_it =
             ids.begin(); id_it != ids.end(); ++id_it)
      {
        removeMatchingItems(id_it->getHits(), decoy_filter);
      }
    }

    
    /**
       @brief Filters peptide or protein identifications according to the given proteins

       Hits with no matching protein accession in @p accessions are removed.

       @note The ranks of the hits may be invalidated.
    */
    template <class IdentificationType>
    static void keepHitsMatchingProteins(std::vector<IdentificationType>& ids,
                                         const std::set<String> accessions)
    {
      struct HasMatchingAccession<typename IdentificationType::HitType>
        acc_filter(accessions);
      for (typename std::vector<IdentificationType>::iterator id_it =
             ids.begin(); id_it != ids.end(); ++id_it)
      {
        keepMatchingItems(id_it->getHits(), acc_filter);
      }
    }


    /// Returns the total number of peptide/protein hits in a vector of peptide/protein identifications
    template <class IdentificationType>
    static Size countHits(const std::vector<IdentificationType>& ids)
    {
      Size counter = 0;
      for (typename std::vector<IdentificationType>::const_iterator id_it =
             ids.begin(); id_it != ids.end(); ++id_it)
      {
        counter += id_it->getHits().size();
      }
      return counter;
    }


    /// Updates the hit ranks on all peptide or protein IDs
    template <class IdentificationType>
    static void updateHitRanks(std::vector<IdentificationType>& ids)
    {
      for (typename std::vector<IdentificationType>::iterator it = ids.begin();
           it != ids.end(); ++it)
      {
        it->assignRanks();
      }
    }


    /**
       @brief Extracts all unique peptide sequences from a list of peptide IDs
       
       @param peptides Input
       @param sequences Output
       @param ignore_mods Extract sequences without modifications?
    */
    static void extractPeptideSequences(
      const std::vector<PeptideIdentification>& peptides, 
      std::set<String>& sequences, bool ignore_mods);


    /**
       @brief Remove duplicate peptide hits from a peptide identification, keeping only unique hits.

       Hits are considered duplicated if they compare as equal using PeptideHit::operator== (i.e. not only the sequences have to match!).
    */
    static void removeDuplicatePeptideHits(std::vector<PeptideIdentification>&
                                           peptides);


    /**
       @brief Filter peptide identifications according to mass deviation.

       Only peptide hits with a low mass deviation (between theoretical peptide mass and precursor mass) are kept.

       @param identification Input/output
       @param mass_error Threshold for the mass deviation
       @param unit_ppm Is @p mass_error given in PPM?

       @note The ranks of the hits may be invalidated.
    */
    static void filterPeptidesByMZError(
      std::vector<PeptideIdentification>& peptides, double mass_error, 
      bool unit_ppm);


    /**
       @brief Filters peptide identifications keeping only the single best-scoring hit per ID.

       @param peptides Input/output
       @param strict If set, keep the best hit only if its score is unique - i.e. ties are not allowed. (Otherwise the first hit with the best score is kept.)
    */
    static void keepBestPeptideHits(
      std::vector<PeptideIdentification>& peptides, bool strict = false);


    /**
       @brief Filters peptide identifications according to peptide sequence length.

       Only peptide hits with a sequence length between @p min_length and @p max_length (both inclusive) are kept.
       @p max_length is ignored if it is smaller than @p min_length.

       @note The ranks of the hits may be invalidated.
    */
    static void filterPeptidesByLength(
      std::vector<PeptideIdentification>& peptides, Size min_length,
      Size max_length = UINT_MAX);


    /**
       @brief Filters peptide identifications according to charge state.

       Only peptide hits with a charge state between @p min_charge and @p max_charge (both inclusive) are kept.
       @p max_charge is ignored if it is smaller than @p min_charge.

       @note The ranks of the hits may be invalidated.
    */
    static void filterPeptidesByCharge(
      std::vector<PeptideIdentification>& peptides, Int min_charge,
      Int max_charge);


    /// Keeps only peptide hits that have at least one of the given modifications
    static void keepPeptidesWithMatchingModifications(
      std::vector<PeptideIdentification>& peptides,
      const std::set<String>& modifications);


    /**
       @brief Removes all peptide hits with a sequence that matches one in @p bad_peptides.

       If @p ignore_mods is set, unmodified sequences are generated and compared to the given ones.

       @note The ranks of the hits may be invalidated.
    */
    static void removePeptidesWithMatchingSequences(
      std::vector<PeptideIdentification>& peptides,
      const std::vector<PeptideIdentification>& bad_peptides, bool ignore_mods);


	  /**
       @brief Filters peptide identifications according to p-values from RTPredict.

       Filters the peptide hits by the probability (p-value) of a correct peptide identification having a deviation between observed and predicted RT equal to or greater than allowed.

       @note The ranks of the hits may be invalidated.
    */
    static void filterPeptidesByRTPredictPValue(
      std::vector<PeptideIdentification>& peptides,
      const String& metavalue_key, double threshold = 0.05);


    /// Removes all peptides that are not annotated as unique for a protein (by PeptideIndexer)
    static void keepUniquePeptidesPerProtein(std::vector<PeptideIdentification>&
                                             peptides);


    /// Removes protein hits from @p proteins that are not referenced by a peptide in @p peptides
    static void removeUnreferencedProteins(
      std::vector<ProteinIdentification>& proteins,
      const std::vector<PeptideIdentification>& peptides);


    /**
       @brief Removes references to missing proteins

       Only PeptideEvidence entries that reference protein hits in @p proteins are kept in the peptide hits.

       If @p remove_peptides_without_reference is set, peptide hits without any remaining protein reference are removed.
    */
    static void updateProteinReferences(
      std::vector<PeptideIdentification>& peptides,
      const std::vector<ProteinIdentification>& proteins,
      bool remove_peptides_without_reference = false);


    /// Filters peptide identifications by precursor RT, keeping only IDs in the given range
    static void filterPeptidesByRT(std::vector<PeptideIdentification>& peptides,
                                   double min_rt, double max_rt);


    /// Filters peptide identifications by precursor m/z, keeping only IDs in the given range
    static void filterPeptidesByMZ(std::vector<PeptideIdentification>& peptides,
                                   double min_mz, double max_mz);


    /**
       @brief Update protein groups after protein hits were filtered

       @param groups Input/output protein groups
       @param hits Available protein hits (all others are removed from the groups)

       @return Returns whether the groups are still valid (which is the case if only whole groups, if any, were removed).
    */
    static bool updateProteinGroups(
      std::vector<ProteinIdentification::ProteinGroup>& groups,
      const std::vector<ProteinHit>& hits);


    /// Filters an MS/MS experiment according to fractions of the significance thresholds
    template <class PeakT>
    static void filterHitsBySignificance(MSExperiment<PeakT>& experiment,
                                         double peptide_threshold_fraction,
                                         double protein_threshold_fraction)
    {
      // filter protein hits:
      filterHitsBySignificance(experiment.getProteinIdentifications(),
                               protein_threshold_fraction);
      // don't remove empty protein IDs - they contain search meta data and may
      // be referenced by peptide IDs (via run ID)

      // filter peptide hits:
      for (typename MSExperiment<PeakT>::Iterator exp_it = experiment.begin();
           exp_it != experiment.end(); ++exp_it)
      {
        filterHitsBySignificance(exp_it->getPeptideIdentifications(),
                                 peptide_threshold_fraction);
        removeEmptyIdentifications(exp_it->getPeptideIdentifications());
        updateProteinReferences(exp_it->getPeptideIdentifications(), 
                                experiment.getProteinIdentifications());
      }
      // @TODO: remove proteins that aren't referenced by peptides any more?
    }


    /// Filters an MS/MS experiment according to score thresholds
    template <class PeakT>
    static void filterHitsByScore(MSExperiment<PeakT>& experiment,
                                  double peptide_threshold_score,
                                  double protein_threshold_score)
    {
      // filter protein hits:
      filterHitsByScore(experiment.getProteinIdentifications(),
                        protein_threshold_score);
      // don't remove empty protein IDs - they contain search meta data and may
      // be referenced by peptide IDs (via run ID)

      // filter peptide hits:
      for (typename MSExperiment<PeakT>::Iterator exp_it = experiment.begin();
           exp_it != experiment.end(); ++exp_it)
      {
        filterHitsByScore(exp_it->getPeptideIdentifications(), 
                          peptide_threshold_score);
        removeEmptyIdentifications(exp_it->getPeptideIdentifications());
        updateProteinReferences(exp_it->getPeptideIdentifications(), 
                                experiment.getProteinIdentifications());
      }
      // @TODO: remove proteins that aren't referenced by peptides any more?
    }


    /// Filters an MS/MS experiment by keeping the N best peptide hits for every spectrum
    template <class PeakT>
    static void keepNBestHits(MSExperiment<PeakT>& experiment, Size n)
    {
      // don't filter the protein hits by "N best" here - filter the peptides
      // and update the protein hits!
      std::vector<PeptideIdentification> all_peptides; // IDs from all spectra

      // filter peptide hits:
      for (typename MSExperiment<PeakT>::Iterator exp_it = experiment.begin();
           exp_it != experiment.end(); ++exp_it)
      {
        std::vector<PeptideIdentification>& peptides = 
          exp_it->getPeptideIdentifications();
        keepNBestHits(peptides, n);
        removeEmptyIdentifications(peptides);
        updateProteinReferences(peptides, 
                                experiment.getProteinIdentifications());
        all_peptides.insert(all_peptides.end(), peptides.begin(),
                            peptides.end());
      }
      // update protein hits:
      removeUnreferencedProteins(experiment.getProteinIdentifications(),
                                 all_peptides);
    }


    /// Filters an MS/MS experiment according to the given proteins
    template <class PeakT>
    static void keepHitsMatchingProteins(
      MSExperiment<PeakT>& experiment, 
      const std::vector<FASTAFile::FASTAEntry>& proteins)
    {
      std::set<String> accessions;
      for (std::vector<FASTAFile::FASTAEntry>::const_iterator it =
             proteins.begin(); it != proteins.end(); ++it)
      {
        accessions.insert(it->identifier);
      }

      // filter protein hits:
      keepHitsMatchingProteins(experiment.getProteinIdentifications(),
                               accessions);
      updateHitRanks(experiment.getProteinIdentifications());

      // filter peptide hits:
      for (typename MSExperiment<PeakT>::Iterator exp_it = experiment.begin();
           exp_it != experiment.end(); ++exp_it)
      {
        if (exp_it->getMSLevel() == 2)
        {
          keepHitsMatchingProteins(exp_it->getPeptideIdentifications(),
                                   accessions);
          removeEmptyIdentifications(exp_it->getPeptideIdentifications());
          updateHitRanks(exp_it->getPeptideIdentifications());
        }
      }
    }

  };

} // namespace OpenMS

#endif // OPENMS_FILTERING_ID_IDFILTER_H
