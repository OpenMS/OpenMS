// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>

#include <algorithm>
#include <climits>
#include <vector>
#include <set>
#include <map>
#include <unordered_set>

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

    /// Typedefs
    typedef std::map<Int, PeptideHit*> ChargeToPepHitP;
    typedef std::unordered_map<std::string, ChargeToPepHitP> SequenceToChargeToPepHitP;
    typedef std::map<std::string, SequenceToChargeToPepHitP> RunToSequenceToChargeToPepHitP;

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

      HasGoodScore(double score_, bool higher_score_better_) :
        score(score_),
        higher_score_better(higher_score_better_)
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

      HasMaxRank(Size rank_):
        rank(rank_)
      {
        if (rank_ == 0)
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

      HasMetaValue(const String& key_, const DataValue& value_):
        key(key_),
        value(value_)
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

      HasMaxMetaValue(const String& key_, const double& value_):
        key(key_),
        value(value_)
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
    struct HasMatchingAccessionUnordered
    {
      typedef HitType argument_type; // for use as a predicate

      const std::unordered_set<String>& accessions;

      HasMatchingAccessionUnordered(const std::unordered_set<String>& accessions_):
        accessions(accessions_)
      {}

      bool operator()(const PeptideHit& hit) const
      {
        for (const auto& it : hit.extractProteinAccessionsSet())
        {
          if (accessions.count(it) > 0) return true;
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
       @brief Given a list of protein accessions, do any occur in the annotation(s) of this hit?

       @note This predicate also works for peptide evidence (class PeptideEvidence).
    */
    template <class HitType>
    struct HasMatchingAccession
    {
      typedef HitType argument_type; // for use as a predicate

      const std::set<String>& accessions;

      HasMatchingAccession(const std::set<String>& accessions_):
        accessions(accessions_)
      {}

      bool operator()(const PeptideHit& hit) const
      {
        for (const auto& it : hit.extractProteinAccessionsSet())
        {
          if (accessions.count(it) > 0) return true;
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
        hits.erase(std::remove_if(hits.begin(), hits.end(), (*this)),
                   hits.end());
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
          OPENMS_LOG_WARN << "Invalid limits! Peptide '" << evidence.getProteinAccession() << "' not filtered" << std::endl;
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
            OPENMS_LOG_WARN << "Peptide accession not available! Skipping Evidence." << std::endl;
          }
          else
          {
            OPENMS_LOG_WARN << "Peptide accession '" << evidence.getProteinAccession()
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

    /// Remove Hit items that satisfy a condition in one of our ID containers (e.g. vector of Peptide or ProteinIDs)
    template <class IDContainer, class Predicate>
    static void removeMatchingItemsUnroll(IDContainer& items, const Predicate& pred)
    {
      for (auto& item : items)
      {
        removeMatchingItems(item.getHits(), pred);
      }
    }

    /// Keep Hit items that satisfy a condition in one of our ID containers (e.g. vector of Peptide or ProteinIDs)
    template <class IDContainer, class Predicate>
    static void keepMatchingItemsUnroll(IDContainer& items, const Predicate& pred)
    {
      for (auto& item : items)
      {
        keepMatchingItems(item.getHits(), pred);
      }
    }

    template <class MapType, class Predicate>
    static void keepMatchingPeptideHits(MapType& prot_and_pep_ids, Predicate& pred)
    {
      for (auto& feat : prot_and_pep_ids)
      {
        keepMatchingItemsUnroll(feat.getPeptideIdentifications(), pred);
      }
      keepMatchingItemsUnroll(prot_and_pep_ids.getUnassignedPeptideIdentifications(), pred);
    }

    template <class MapType, class Predicate>
    static void removeMatchingPeptideHits(MapType& prot_and_pep_ids, Predicate& pred)
    {
      for (auto& feat : prot_and_pep_ids)
      {
        removeMatchingItemsUnroll(feat.getPeptideIdentifications(), pred);
      }
      removeMatchingItemsUnroll(prot_and_pep_ids.getUnassignedPeptideIdentifications(), pred);
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
       @brief Removes references to missing proteins

       Only PeptideEvidence entries that reference protein hits in @p proteins are kept in the peptide hits.

       If @p remove_peptides_without_reference is set, peptide hits without any remaining protein reference are removed.
    */
    static void updateProteinReferences(
        ConsensusMap& cmap,
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
      @brief Filters peptide or protein identifications according to the score of the hits.

      Only peptide/protein hits with a score at least as good as @p threshold_score are kept. Score orientation (are higher scores better?) is taken into account.
    */
    template <class IdentificationType>
    static void filterHitsByScore(IdentificationType& id,
                                  double threshold_score)
    {
        struct HasGoodScore<typename IdentificationType::HitType> score_filter(
            threshold_score, id->isHigherScoreBetter());
        keepMatchingItems(id->getHits(), score_filter);
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
      struct HasMatchingAccession<typename IdentificationType::HitType> acc_filter(accessions);
      for (auto& id_it : ids)
      {
        removeMatchingItems(id_it.getHits(), acc_filter);
      }
    }

    /**
       @brief Filters peptide or protein identifications according to the given proteins (positive).

       Hits with no matching protein accession in @p accessions are removed.

       @note The ranks of the hits may be invalidated.
    */
    template <class IdentificationType>
    static void keepHitsMatchingProteins(std::vector<IdentificationType>& ids,
                                         const std::set<String>& accessions)
    {
      struct HasMatchingAccession<typename IdentificationType::HitType> acc_filter(accessions);
      for (auto& id_it : ids)
      {
        keepMatchingItems(id_it.getHits(), acc_filter);
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

    /// Filters a Consensus/FeatureMap by keeping the N best peptide hits for every spectrum
    template <class MapType>
    static void keepNBestPeptideHits(MapType& map, Size n)
    {
      // The rank predicate needs annotated ranks, not sure if they are always updated. Use the following instead,
      // which sorts Hits first.
      for (auto& feat : map)
      {
        keepNBestHits(feat.getPeptideIdentifications(), n);
      }
      keepNBestHits(map.getUnassignedPeptideIdentifications(), n);
    }

    template <class MapType>
    static void removeEmptyIdentifications(MapType& prot_and_pep_ids)
    {
      removeMatchingPeptideHits(prot_and_pep_ids, HasNoHits<PeptideHit>());
    }

    /// Filters PeptideHits from PeptideIdentification by keeping only the best peptide hits for every peptide sequence
    static void keepBestPerPeptide(std::vector<PeptideIdentification>& pep_ids, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
    {
      annotateBestPerPeptide(pep_ids, ignore_mods, ignore_charges, nr_best_spectrum);
      HasMetaValue<PeptideHit> best_per_peptide{"best_per_peptide", 1};
      keepMatchingItemsUnroll(pep_ids, best_per_peptide);
    }

    static void keepBestPerPeptidePerRun(std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
    {
      annotateBestPerPeptidePerRun(prot_ids, pep_ids, ignore_mods, ignore_charges, nr_best_spectrum);
      HasMetaValue<PeptideHit> best_per_peptide{"best_per_peptide", 1};
      keepMatchingItemsUnroll(pep_ids, best_per_peptide);
    }

    template <class MapType>
    static void keepBestPerPeptidePerRun(MapType& prot_and_pep_ids, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
    {
      const auto& prot_ids = prot_and_pep_ids.getProteinIdentifications();

      RunToSequenceToChargeToPepHitP best_peps_per_run;
      for (const auto& idrun : prot_ids)
      {
        best_peps_per_run[idrun.getIdentifier()] = SequenceToChargeToPepHitP();
      }

      for (auto& feat : prot_and_pep_ids)
      {
        annotateBestPerPeptidePerRunWithData(best_peps_per_run, feat.getPeptideIdentifications(), ignore_mods, ignore_charges, nr_best_spectrum);
      }

      annotateBestPerPeptidePerRunWithData(best_peps_per_run, prot_and_pep_ids.getUnassignedPeptideIdentifications(), ignore_mods, ignore_charges, nr_best_spectrum);

      HasMetaValue<PeptideHit> best_per_peptide{"best_per_peptide", 1};
      keepMatchingPeptideHits(prot_and_pep_ids, best_per_peptide);
    }

    /// Annotates PeptideHits from PeptideIdentification if it is the best peptide hit for its peptide sequence
    /// Adds metavalue "bestForItsPeps" which can be used for additional filtering.
    static void annotateBestPerPeptidePerRun(const std::vector<ProteinIdentification>& prot_ids, std::vector<PeptideIdentification>& pep_ids, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
    {
      RunToSequenceToChargeToPepHitP best_peps_per_run;
      for (const auto& id : prot_ids)
      {
        best_peps_per_run[id.getIdentifier()] = SequenceToChargeToPepHitP();
      }
      annotateBestPerPeptidePerRunWithData(best_peps_per_run, pep_ids, ignore_mods, ignore_charges, nr_best_spectrum);
    }

    /// Annotates PeptideHits from PeptideIdentification if it is the best peptide hit for its peptide sequence
    /// Adds metavalue "bestForItsPeps" which can be used for additional filtering.
    /// To be used when a RunToSequenceToChargeToPepHitP map is already available
    static void annotateBestPerPeptidePerRunWithData(RunToSequenceToChargeToPepHitP& best_peps_per_run, std::vector<PeptideIdentification>& pep_ids, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
    {
      for (auto &pep : pep_ids)
      {
        SequenceToChargeToPepHitP& best_pep = best_peps_per_run[pep.getIdentifier()];
        annotateBestPerPeptideWithData(best_pep, pep, ignore_mods, ignore_charges, nr_best_spectrum);
      }
    }

    /// Annotates PeptideHits from PeptideIdentification if it is the best peptide hit for its peptide sequence
    /// Adds metavalue "bestForItsPeps" which can be used for additional filtering.
    /// Does not check Run information and just goes over all Peptide IDs
    static void annotateBestPerPeptide(std::vector<PeptideIdentification>& pep_ids, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
    {
      SequenceToChargeToPepHitP best_pep;
      for (auto& pep : pep_ids)
      {
        annotateBestPerPeptideWithData(best_pep, pep, ignore_mods, ignore_charges, nr_best_spectrum);
      }
    }

    /// Annotates PeptideHits from PeptideIdentification if it is the best peptide hit for its peptide sequence
    /// Adds metavalue "bestForItsPeps" which can be used for additional filtering.
    /// Does not check Run information and just goes over all Peptide IDs
    /// To be used when a SequenceToChargeToPepHitP map is already available
    static void annotateBestPerPeptideWithData(SequenceToChargeToPepHitP& best_pep, PeptideIdentification& pep, bool ignore_mods, bool ignore_charges, Size nr_best_spectrum)
    {
        bool higher_score_better = pep.isHigherScoreBetter();
        //make sure that first = best hit
        pep.sort();

        auto pepIt = pep.getHits().begin();
        auto pepItEnd = nr_best_spectrum == 0 || pep.getHits().size() <= nr_best_spectrum ? pep.getHits().end() : pep.getHits().begin() + nr_best_spectrum;
        for (; pepIt != pepItEnd; ++pepIt)
        {
          PeptideHit &hit = *pepIt;

          String lookup_seq;
          if (ignore_mods)
          {
            lookup_seq = hit.getSequence().toUnmodifiedString();
          }
          else
          {
            lookup_seq = hit.getSequence().toString();
          }

          int lookup_charge = 0;
          if (!ignore_charges)
          {
            lookup_charge = hit.getCharge();
          }

          // try to insert
          auto it_inserted = best_pep.emplace(std::move(lookup_seq), ChargeToPepHitP());
          auto it_inserted_chg = it_inserted.first->second.emplace(lookup_charge, &hit);

          PeptideHit* &p = it_inserted_chg.first->second; //now this gets either the old one if already present, or this
          if (!it_inserted_chg.second) //was already present -> possibly update
          {
            if (
                (higher_score_better && (hit.getScore() > p->getScore())) ||
                (!higher_score_better && (hit.getScore() < p->getScore()))
                )
            {
              p->setMetaValue("best_per_peptide", 0);
              hit.setMetaValue("best_per_peptide", 1);
              p = &hit;
            }
            else //note that this was def. not the best
            {
              // TODO if it is only about filtering, we can omit writing this metavalue (absence = false)
              hit.setMetaValue("best_per_peptide", 0);
            }
          }
          else //newly inserted -> first for that sequence (and optionally charge)
          {
            hit.setMetaValue("best_per_peptide", 1);
          }
        }
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


    /// @name Filter functions for class IdentificationData
    ///@{
    static void keepBestMatchPerQuery(
      IdentificationData& id_data,
      IdentificationData::ScoreTypeRef score_ref);

    static void filterQueryMatchesByScore(
      IdentificationData& id_data,
      IdentificationData::ScoreTypeRef score_ref, double cutoff);

    static void removeDecoys(IdentificationData& id_data);
    ///@}

  };

} // namespace OpenMS

