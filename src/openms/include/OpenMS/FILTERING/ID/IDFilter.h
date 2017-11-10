// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <algorithm>
#include <climits>
#include <vector>
#include <set>
#include <map>

namespace OpenMS
{
  /**
    @brief Collection of functions for filtering peptide and protein identifications.

    This class provides functions for filtering collections of peptide or protein identifications according to various criteria.
    It also contains helper functions and classes (functors that implement predicates) that are used in this context.

    The filter functions modify their inputs, rather than creating filtered copies.

    Most filters work on the hit level, i.e. they remove peptide or protein hits from peptide or protein identifications (IDs).
    A few filters work on the ID level instead, i.e. they remove peptide or protein IDs from vectors thereof.
    Independent of this, the inputs for all filter functions are vectors of IDs, because the data most often comes in this form.
    This design also allows many helper objects to be set up only once per vector, rather than once per ID.

    The filter functions for vectors of peptide/protein IDs do not include clean-up steps (e.g. removal of IDs without hits, reassignment of hit ranks, ...).
    They only carry out their specific filtering operations.
    This is so filters can be chained without having to repeat clean-up operations.
    The group of clean-up functions provides helpers that are useful to ensure data integrity after filters have been applied, but it is up to the individual developer to use them when necessary.

    The filter functions for MS/MS experiments do include clean-up steps, because they filter peptide and protein IDs in conjunction and potential contradictions between the two must be eliminated.
  */
  class OPENMS_DLLAPI IDFilter
  {
public:

    /// Constructor
    IDFilter();

    /// Destructor
    virtual ~IDFilter();


    /**
       @name Predicates for peptide or protein hits

       These functors test for some property of a peptide or protein hit
    */
    ///@{

    /// Is the score of this hit at least as good as the given value?
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
       @brief Is the rank of this hit below or at the given cut-off?

       Ranks are counted from one (best), so zero is not a valid cut-off.
    */
    template <class HitType>
    struct HasMaxRank
    {
      typedef HitType argument_type; // for use as a predicate

      Size rank;

      HasMaxRank(Size rank):
        rank(rank)
      {
        if (rank == 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The cut-off value for rank filtering must not be zero!");
        }
      }

      bool operator()(const HitType& hit) const
      {
        Size hit_rank = hit.getRank();
        if (hit_rank == 0)
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No rank assigned to peptide or protein hit");
        }
        return hit_rank <= rank;
      }
    };

    /**
       @brief Is a meta value with given key and value set on this hit?

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

    /// Does a meta value of this hit have at most the given value?
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

    /// Is this a decoy hit?
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

    /**
       @brief Given a list of protein accessions, do any occur in the annotation(s) of this hit?

       @note This predicate also works for peptide evidence (class PeptideEvidence).
    */
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
        std::set<String> present_accessions = hit.extractProteinAccessionsSet();
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

    /**
       @brief Builds a map index of data that have a String index to find matches and return the objects

       @note Currently implemented for FastaEntries and Peptide Evidences
    */
    template <class HitType, class Entry>
    struct GetMatchingItems
    {
      typedef HitType argument_type; // for use as a predicate
      typedef std::map<String, Entry*> ItemMap;//Store pointers to avoid copying data
      ItemMap items;

      GetMatchingItems(std::vector<Entry>& records)
      {
        for(typename std::vector<Entry>::iterator rec_it = records.begin();
            rec_it != records.end(); ++rec_it)
        {
          items[getKey(*rec_it)] = &(*rec_it);
        }
      }

      GetMatchingItems(){}

      const String& getKey(const FASTAFile::FASTAEntry& entry) const
      {
        return entry.identifier;
      }

      bool exists(const HitType& hit) const
      {
        return items.count(getHitKey(hit)) > 0;
      }

      const String& getHitKey(const PeptideEvidence& p) const
      {
        return p.getProteinAccession();
      }

      const Entry& getValue(const PeptideEvidence& evidence) const
      {
        if(!exists(evidence)){
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Accession: '"+ getHitKey(evidence) + "'. peptide evidence accession not in data");
        }
        return *(items.find(getHitKey(evidence))->second);
      }

    };


    ///@}


    /**
       @name Predicates for peptide hits only

       These functors test for some property of peptide hits
    */
    ///@{

    /// Does the sequence of this peptide hit have at least the given length?
    struct HasMinPeptideLength;

    /// Does the charge of this peptide hit have at least the given value?
    struct HasMinCharge;

    /// Is the m/z error of this peptide hit below the given value?
    struct HasLowMZError;

    /**
       @brief Given a list of modifications, do any occur in the sequence of this peptide hit?

       If the list of modifications is empty, return true if the sequence is modified at all.
    */
    struct HasMatchingModification;

    /**
       @brief Is the sequence of this peptide hit among a list of given sequences?

       With @p ignore_mods, the sequence without modifications is compared.
    */
    struct HasMatchingSequence;

    /// Is the list of peptide evidences of this peptide hit empty?
    struct HasNoEvidence;

    
    /**
       @brief Filter Peptide Hit by its digestion product

    */

    class PeptideDigestionFilter
    {
     private:
      EnzymaticDigestion& digestion_;
      Int min_cleavages_;
      Int max_cleavages_;

     public:
      typedef PeptideHit argument_type;
      PeptideDigestionFilter(EnzymaticDigestion& digestion, Int min, Int max) :
      digestion_(digestion), min_cleavages_(min), max_cleavages_(max)
      {}

      static inline Int disabledValue(){ return -1; }
      
      /// Filter function on min max cutoff values to be used with remove_if 
      /// returns true if peptide should be removed (does not pass filter)
      bool operator()(PeptideHit& p)
      {
        return digestion_.filterByMissedCleavages(
          p.getSequence().toUnmodifiedString(),
          [&](const Int missed_cleavages)
          {

            bool max_filter = max_cleavages_ != disabledValue() ? 
                              missed_cleavages > max_cleavages_ : false;
            bool min_filter = min_cleavages_ != disabledValue() ?
                              missed_cleavages < min_cleavages_ : false;
            return max_filter || min_filter;
          });
      }

      void filterPeptideSequences(std::vector<PeptideHit>& hits)
      {
        hits.erase(std::remove_if(hits.begin(), hits.end(), (*this)), hits.end());
      }

    };


    /**
       @brief Is peptide evidence digestion product of some protein

       Keeps all valid products
     */
    struct DigestionFilter
    {
      typedef PeptideEvidence argument_type;

      // Build an accession index to avoid the linear search cost
      GetMatchingItems<PeptideEvidence, FASTAFile::FASTAEntry>accession_resolver_;
      ProteaseDigestion& digestion_;
      bool ignore_missed_cleavages_;
      bool methionine_cleavage_;

      DigestionFilter(std::vector<FASTAFile::FASTAEntry>& entries,
                      ProteaseDigestion& digestion,
                      bool ignore_missed_cleavages,
                      bool methionine_cleavage) :
        accession_resolver_(entries),
        digestion_(digestion),
        ignore_missed_cleavages_(ignore_missed_cleavages),
        methionine_cleavage_(methionine_cleavage)
      {}

      bool operator()(const PeptideEvidence& evidence) const
      {
        if(!evidence.hasValidLimits())
        {
          LOG_WARN << "Invalid limits! Peptide '" << evidence.getProteinAccession() << "' not filtered" << std::endl;
          return true;
        }

        if (accession_resolver_.exists(evidence))
        {
          return digestion_.isValidProduct(
            AASequence::fromString(accession_resolver_.getValue(evidence).sequence),
            evidence.getStart(), evidence.getEnd() - evidence.getStart(), ignore_missed_cleavages_, methionine_cleavage_);
        }
        else
        {
          if (evidence.getProteinAccession().empty())
          {
            LOG_WARN << "Peptide accession not available! Skipping Evidence." << std::endl;
          }
          else
          {
            LOG_WARN << "Peptide accession '" << evidence.getProteinAccession()
                     << "' not found in fasta file!" << std::endl;
          }
          return true;
        }
      }

      void filterPeptideEvidences(std::vector<PeptideIdentification>& peptides)
      {
        IDFilter::FilterPeptideEvidences<IDFilter::DigestionFilter>(*this,peptides);
      }

    };


    ///@}


    /// @name Predicates for peptide or protein identifications
    ///@{

    /// Is the list of hits of this peptide/protein ID empty?
    template <class IdentificationType>
    struct HasNoHits
    {
      typedef IdentificationType argument_type; // for use as a predicate

      bool operator()(const IdentificationType& id) const
      {
        return id.getHits().empty();
      }
    };

    ///@}


    /// @name Predicates for peptide identifications only
    ///@{

    /// Is the retention time of this peptide ID in the given range?
    struct HasRTInRange;

    /// Is the precursor m/z value of this peptide ID in the given range?
    struct HasMZInRange;

    ///@}


    /**
       @name Higher-order filter functions

       Functions for filtering a container based on a predicate
    */
    ///@{

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

    ///@}


    /// @name Helper functions
    ///@{

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

    /**
       @brief Finds the best-scoring hit in a vector of peptide or protein identifications.

       If there are several hits with the best score, the first one is taken.

       @param identifications Vector of peptide or protein IDs, each containing one or more (peptide/protein) hits
       @param assume_sorted Are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at

       @except Exception::InvalidValue if the IDs have different score types (i.e. scores cannot be compared)

       @return true if a hit was present, false otherwise
    */
    template <class IdentificationType>
    static bool getBestHit(
      const std::vector<IdentificationType>& identifications,
      bool assume_sorted, typename IdentificationType::HitType& best_hit)
    {
      if (identifications.empty()) return false;

      typename std::vector<IdentificationType>::const_iterator best_id_it =
        identifications.end();
      typename std::vector<typename IdentificationType::HitType>::const_iterator
        best_hit_it;

      for (typename std::vector<IdentificationType>::const_iterator id_it =
             identifications.begin(); id_it != identifications.end(); ++id_it)
      {
        if (id_it->getHits().empty()) continue;

        if (best_id_it == identifications.end()) // no previous "best" hit
        {
          best_id_it = id_it;
          best_hit_it = id_it->getHits().begin();
        }
        else if (best_id_it->getScoreType() != id_it->getScoreType())
        {
          throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Can't compare scores of different types", best_id_it->getScoreType() + "/" + id_it->getScoreType());
        }

        bool higher_better = best_id_it->isHigherScoreBetter();
        for (typename std::vector<typename IdentificationType::HitType>::
               const_iterator hit_it = id_it->getHits().begin(); hit_it !=
               id_it->getHits().end(); ++hit_it)
        {
          if ((higher_better && (hit_it->getScore() >
                                 best_hit_it->getScore())) ||
              (!higher_better && (hit_it->getScore() <
                                  best_hit_it->getScore())))
          {
            best_hit_it = hit_it;
          }
          if (assume_sorted) break; // only consider the first hit
        }
      }

      if (best_id_it == identifications.end())
      {
        return false; // no hits in any IDs
      }

      best_hit = *best_hit_it;
      return true;
    }

    /**
       @brief Extracts all unique peptide sequences from a list of peptide IDs

       @param peptides Input
       @param sequences Output
       @param ignore_mods Extract sequences without modifications?
    */
    static void extractPeptideSequences(
      const std::vector<PeptideIdentification>& peptides,
      std::set<String>& sequences, bool ignore_mods = false);

    /**
       @brief remove peptide evidences based on a filter
       
       @param filter filter function that overloads ()(PeptideEvidence&) operator
       @param peptides a collection of peptide evidences
     */
    template<class EvidenceFilter>
    static void FilterPeptideEvidences(
      EvidenceFilter& filter,
      std::vector<PeptideIdentification>& peptides)
    {
      for(std::vector<PeptideIdentification>::iterator pep_it = peptides.begin();
          pep_it != peptides.end(); ++pep_it)
      {
        for(std::vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
            hit_it != pep_it->getHits().end(); ++hit_it )
        {
          std::vector<PeptideEvidence> evidences;
          remove_copy_if(hit_it->getPeptideEvidences().begin(),
                         hit_it->getPeptideEvidences().end(),
                         back_inserter(evidences), 
                         std::not1(filter));
          hit_it->setPeptideEvidences(evidences);
        }
      }
    }
    

    ///@}


    /// @name Clean-up functions
    ///@{

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

    /**
       @brief Update protein groups after protein hits were filtered

       @param groups Input/output protein groups
       @param hits Available protein hits (all others are removed from the groups)

       @return Returns whether the groups are still valid (which is the case if only whole groups, if any, were removed).
    */
    static bool updateProteinGroups(
      std::vector<ProteinIdentification::ProteinGroup>& groups,
      const std::vector<ProteinHit>& hits);

    ///@}


    /// @name Filter functions for peptide or protein IDs
    ///@{

    /// Removes peptide or protein identifications that have no hits in them
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
       @brief Filters peptide or protein identifications according to the significance threshold of the hits.

       Only peptide/protein hits which reach a score above (or below, depending on score orientation) @p threshold_fraction * @p significance_threshold (as stored in the ID) are kept.
    */
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

       The hits between @p min_rank and @p max_rank (both inclusive) in each ID are kept.
       Counting starts at 1, i.e. the best (highest/lowest scoring) hit has rank 1.
       The ranks are (re-)computed before filtering.
       @p max_rank is ignored if it is smaller than @p min_rank.

       Note that there may be several hits with the same rank in a peptide or protein ID (if the scores are the same).

       This method is useful if a range of higher hits is needed for decoy fairness analysis.

       @note The ranks of the hits may be invalidated.
    */
    template <class IdentificationType>
    static void filterHitsByRank(std::vector<IdentificationType>& ids,
                                 Size min_rank, Size max_rank)
    {
      updateHitRanks(ids);
      if (min_rank > 1)
      {
        struct HasMaxRank<typename IdentificationType::HitType>
          rank_filter(min_rank - 1);
        for (typename std::vector<IdentificationType>::iterator id_it =
               ids.begin(); id_it != ids.end(); ++id_it)
        {
          removeMatchingItems(id_it->getHits(), rank_filter);
        }
      }
      if (max_rank >= min_rank)
      {
        struct HasMaxRank<typename IdentificationType::HitType>
          rank_filter(max_rank);
        for (typename std::vector<IdentificationType>::iterator id_it =
               ids.begin(); id_it != ids.end(); ++id_it)
        {
          keepMatchingItems(id_it->getHits(), rank_filter);
        }
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
       @brief Filters peptide or protein identifications according to the given proteins (negative).

       Hits with a matching protein accession in @p accessions are removed.

       @note The ranks of the hits may be invalidated.
    */
    template <class IdentificationType>
    static void removeHitsMatchingProteins(std::vector<IdentificationType>& ids,
                                           const std::set<String> accessions)
    {
      struct HasMatchingAccession<typename IdentificationType::HitType>
        acc_filter(accessions);
      for (typename std::vector<IdentificationType>::iterator id_it =
             ids.begin(); id_it != ids.end(); ++id_it)
      {
        removeMatchingItems(id_it->getHits(), acc_filter);
      }
    }

    /**
       @brief Filters peptide or protein identifications according to the given proteins (positive).

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

    

    ///@}


    /// @name Filter functions for peptide IDs only
    ///@{

    /**
       @brief Filters peptide identifications keeping only the single best-scoring hit per ID.

       @param peptides Input/output
       @param strict If set, keep the best hit only if its score is unique - i.e. ties are not allowed. (Otherwise all hits with the best score is kept.)
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

    /// Filters peptide identifications by precursor RT, keeping only IDs in the given range
    static void filterPeptidesByRT(std::vector<PeptideIdentification>& peptides,
                                   double min_rt, double max_rt);

    /// Filters peptide identifications by precursor m/z, keeping only IDs in the given range
    static void filterPeptidesByMZ(std::vector<PeptideIdentification>& peptides,
                                   double min_mz, double max_mz);

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
       @brief Digest a collection of proteins and filter PeptideEvidences based on specificity
       PeptideEvidences of peptides are removed if the digest of a protein did not produce the peptide sequence
       @param filter filter function on PeptideEvidence level
       @param peptides PeptideIdentification that will be scanned and filtered
     */
    template <class Filter>
    static void filterPeptideEvidences(
      Filter& filter,
      std::vector<PeptideIdentification>& peptides);

	  /**
       @brief Filters peptide identifications according to p-values from RTPredict.

       Filters the peptide hits by the probability (p-value) of a correct peptide identification having a deviation between observed and predicted RT equal to or greater than allowed.

       @param peptides Input/output
       @param metavalue_key Name of the meta value that holds the p-value: "predicted_RT_p_value" or "predicted_RT_p_value_first_dim"
       @param threshold P-value threshold

       @note The ranks of the hits may be invalidated.
    */
    static void filterPeptidesByRTPredictPValue(
      std::vector<PeptideIdentification>& peptides,
      const String& metavalue_key, double threshold = 0.05);

    /// Removes all peptide hits that have at least one of the given modifications
    static void removePeptidesWithMatchingModifications(
      std::vector<PeptideIdentification>& peptides,
      const std::set<String>& modifications);

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
      const std::vector<PeptideIdentification>& bad_peptides,
      bool ignore_mods = false);

    /**
       @brief Removes all peptide hits with a sequence that does not match one in @p good_peptides.

       If @p ignore_mods is set, unmodified sequences are generated and compared to the given ones.

       @note The ranks of the hits may be invalidated.
    */
    static void keepPeptidesWithMatchingSequences(
      std::vector<PeptideIdentification>& peptides,
      const std::vector<PeptideIdentification>& good_peptides,
      bool ignore_mods = false);

   /// Removes all peptides that are not annotated as unique for a protein (by PeptideIndexer)
    static void keepUniquePeptidesPerProtein(std::vector<PeptideIdentification>&
                                             peptides);

    /**
       @brief Removes duplicate peptide hits from each peptide identification, keeping only unique hits (per ID).

       By default, hits are considered duplicated if they compare as equal using PeptideHit::operator==. However, if @p seq_only is set, only the sequences (incl. modifications) are compared. In both cases, the first occurrence of each hit in a peptide ID is kept, later ones are removed.
    */
    static void removeDuplicatePeptideHits(std::vector<PeptideIdentification>&
                                           peptides, bool seq_only = false);

 ///@}


    /// @name Filter functions for MS/MS experiments
    ///@{

    /// Filters an MS/MS experiment according to score thresholds
    static void filterHitsByScore(PeakMap& experiment,
                                  double peptide_threshold_score,
                                  double protein_threshold_score)
    {
      // filter protein hits:
      filterHitsByScore(experiment.getProteinIdentifications(),
                        protein_threshold_score);
      // don't remove empty protein IDs - they contain search meta data and may
      // be referenced by peptide IDs (via run ID)

      // filter peptide hits:
      for (PeakMap::Iterator exp_it = experiment.begin();
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

    /// Filters an MS/MS experiment according to fractions of the significance thresholds
    static void filterHitsBySignificance(PeakMap& experiment,
                                         double peptide_threshold_fraction,
                                         double protein_threshold_fraction)
    {
      // filter protein hits:
      filterHitsBySignificance(experiment.getProteinIdentifications(),
                               protein_threshold_fraction);
      // don't remove empty protein IDs - they contain search meta data and may
      // be referenced by peptide IDs (via run ID)

      // filter peptide hits:
      for (PeakMap::Iterator exp_it = experiment.begin();
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

    /// Filters an MS/MS experiment by keeping the N best peptide hits for every spectrum
    static void keepNBestHits(PeakMap& experiment, Size n)
    {
      // don't filter the protein hits by "N best" here - filter the peptides
      // and update the protein hits!
      std::vector<PeptideIdentification> all_peptides; // IDs from all spectra

      // filter peptide hits:
      for (PeakMap::Iterator exp_it = experiment.begin();
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
    static void keepHitsMatchingProteins(
      PeakMap& experiment,
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
      for (PeakMap::Iterator exp_it = experiment.begin();
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

    ///@}


  };

} // namespace OpenMS

#endif // OPENMS_FILTERING_ID_IDFILTER_H
