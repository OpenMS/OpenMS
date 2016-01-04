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


    /// gets the best scoring peptide hit from a vector of peptide identifications
    /// @param identifications Vector of peptide ids, each containing one or more peptide hits
    /// @param assume_sorted are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
    /// @return true if a hit was present, false otherwise
    template <class IdentificationType>
    static bool getBestHit(const std::vector<IdentificationType> identifications, bool assume_sorted, PeptideHit& best_hit)
    {
      if (identifications.size() == 0) return false;

      bool is_higher_score_better = identifications[0].isHigherScoreBetter();
      double best_score = (is_higher_score_better ? -1 : 1) * std::numeric_limits<double>::max(); // worst score we can think of

      Size best_i_index(0), best_h_index(0);
      Size max_h(-1);
      bool hit_found(false);
      // determine best scoring hit
      for (Size i = 0; i != identifications.size(); ++i)
      {
        if (identifications[i].getHits().size() == 0) continue; // empty hits
        
        is_higher_score_better = identifications[i].isHigherScoreBetter();
        max_h = (assume_sorted ? 1 : identifications[i].getHits().size());
        hit_found = true;
        for (Size h = 0; h < max_h; ++h)
        {
          double score = identifications[i].getHits()[h].getScore();
          // better score?
          if (score > best_score * (is_higher_score_better ? 1 : -1))
          {
            best_score = score;
            best_i_index = i;
            best_h_index = h;
          }
        }
      }

      if (!hit_found)
      {
        return false;// all hits were empty 
      }

      best_hit = identifications[best_i_index].getHits()[best_h_index];
      return true;

    }


    /// filters a ProteinIdentification or PeptideIdentification by only allowing peptides/proteins which reach a score above (or below, depending on score orientation) @p threshold_fraction * SignificanceThreshold
    template <class IdentificationType>
    static void filterIdentificationsByThreshold(
      const IdentificationType& identification, double threshold_fraction,
      IdentificationType& filtered_identification)
    {
      double threshold_score = (threshold_fraction *
                                identification.getSignificanceThreshold());
      
      filterIdentificationsByScore<IdentificationType>(
        identification, threshold_score, filtered_identification);
    }


    /**
      @brief filters a ProteinIdentification or PeptideIdentification corresponding to the @p threshold_score

      If the method higherScoreBetter() returns true for the IdentificationType all hits with a score
      smaller than @p threshold_score are removed. Otherwise all hits with a score bigger than
      @p threshold_score are removed.
    */
    template <class IdentificationType>
    static void filterIdentificationsByScore(const IdentificationType& identification, double threshold_score, IdentificationType& filtered_identification)
    {
      filtered_identification = identification;
      
      bool higher_better = filtered_identification.isHigherScoreBetter();
      struct HasGoodScore<typename IdentificationType::HitType> score_filter(
        threshold_score, higher_better);
      keepMatchingItems(filtered_identification.getHits(), score_filter);

      filtered_identification.assignRanks();
    }


    /**
      @brief filters a ProteinIdentification or PeptideIdentification corresponding to the score.

      If the method higherScoreBetter() returns true for the IdentificationType the
      @p n highest scoring hits are kept. Otherwise the @p n lowest scoring hits are kept.
    */
    template <class IdentificationType>
    static void filterIdentificationsByBestNHits(const IdentificationType& identification, Size n, IdentificationType& filtered_identification)
    {
      filtered_identification = identification;
      filtered_identification.sort();
      if (filtered_identification.getHits().size() > n)
      {
        filtered_identification.getHits().resize(n);
      }

      filtered_identification.assignRanks();
    }


    /**
      @brief filters a ProteinIdentification or PeptideIdentification corresponding to the score.

      If the method higherScoreBetter() returns true for the IdentificationType the
      @p n to @p m highest scoring hits are kept. Otherwise the @p n to @p m lowest scoring hits are kept. Counting starts at 1, i.e. the best (highest/lowest scoring) hit is hit 1; the boundaries @p n and @p m are included in the range.
      This method is useful if a range of higher hits are used for decoy fairness analysis.
    */
    template <class IdentificationType>
    static void filterIdentificationsByBestNToMHits(
      const IdentificationType& identification, Size n, Size m,
      IdentificationType& filtered_identification)
    {
      if (n > m) std::swap(n, m);

      filtered_identification = identification;
      filtered_identification.sort();

      std::vector<PeptideHit>& hits = filtered_identification.getHits();
      Size size = hits.size();
      if (m < size) hits.erase(hits.begin() + m, hits.end());
      if (n > 1) hits.erase(hits.begin(), hits.begin() + n - 1);

      filtered_identification.assignRanks();
    }
    
    
    /**
     @brief filters a ProteinIdentification or PeptideIdentification corresponding to their decoy information.
     
     Checks for "target_decoy" or "isDecoy" meta values and removes protein/peptide hits if the values are "decoy" or "true", respectively.
     */
    template <class IdentificationType>
    static void filterIdentificationsByDecoy(const IdentificationType& identification, IdentificationType& filtered_identification)
    {
      filtered_identification = identification;
      struct HasDecoyAnnotation<typename IdentificationType::HitType> 
        decoy_filter;
      removeMatchingItems(filtered_identification.getHits(), decoy_filter);

      filtered_identification.assignRanks();
    }

    
    /**
       @brief Filters a peptide or protein identification corresponding to the given proteins

       Hits with no matching protein accession in @p fasta are removed.
    */
    template <class IdentificationType>
    static void filterIdentificationsByProteins(
      const IdentificationType& identification, 
      const std::vector<FASTAFile::FASTAEntry>& fasta, 
      IdentificationType& filtered_identification)
    {
      filtered_identification = identification;
      
      std::set<String> accessions;
      for (std::vector<FASTAFile::FASTAEntry>::const_iterator it = 
             fasta.begin(); it != fasta.end(); ++it)
      {
        accessions.insert(it->identifier);
      }
      
      struct HasMatchingAccession<typename IdentificationType::HitType>
        acc_filter(accessions);
      keepMatchingItems(filtered_identification.getHits(), acc_filter);
      
      filtered_identification.assignRanks();
    }


    /**
       @brief Filters a peptide or protein identification corresponding to the given proteins.

       Hits with no matching protein accession in @p proteins are removed.
    */
    template <class IdentificationType>
    static void filterIdentificationsByProteinAccessions(
      const IdentificationType& identification, const StringList& proteins,
      IdentificationType& filtered_identification)
    {
      filtered_identification = identification;
      
      std::set<String> accessions(proteins.begin(), proteins.end());
      struct HasMatchingAccession<typename IdentificationType::HitType>
        acc_filter(accessions);
      keepMatchingItems(filtered_identification.getHits(), acc_filter);
      
      filtered_identification.assignRanks();
    }


    /// filters a PeptideIdentification keeping only the best scoring hits (if @p strict is set, keeping only the best hit only if it is the only hit with that score)
    static void filterIdentificationsByBestHits(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, bool strict = false);

    /// removes all peptide hits having a sequence equal to a String in @p peptides. If @p ignore_modifications is set, the unmodified versions are generated and compared to the set of Strings.
    static void filterIdentificationsByExclusionPeptides(const PeptideIdentification& identification, const std::set<String>& peptides, bool ignore_modifications, PeptideIdentification& filtered_identification);

    /// Only peptides having a length l with @p min_length <= l <= @p max_length will be kept.
    /// @p max_length will be ignored if it is smaller than @p min_length.
    static void filterIdentificationsByLength(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, Size min_length, Size max_length = UINT_MAX);

    /// only peptides that have a charge equal to or greater than @p charge will be kept
    static void filterIdentificationsByCharge(const PeptideIdentification& identification, Int charge, PeptideIdentification& filtered_identification);

    /// only peptides having a variable modification will be kept
    static void filterIdentificationsByVariableModifications(const PeptideIdentification& identification, std::vector<String>& fixed_modifications, PeptideIdentification& filtered_identification);

    /// only protein hits in @p identification which are referenced by a peptide in @p peptide_identifications are kept
    static void removeUnreferencedProteinHits(const ProteinIdentification& identification, const std::vector<PeptideIdentification>& peptide_identifications, ProteinIdentification& filtered_identification);

    /// only peptide hits in @p peptide_identifications which are referenced by a protein in @p identification are kept
    static void removeUnreferencedPeptideHits(const ProteinIdentification& identification, std::vector<PeptideIdentification>& peptide_identifications, bool delete_unreferenced_peptide_hits = false);

    /// if a peptide hit occurs more than once per peptide ID, only one instance is kept
    static void filterIdentificationsUnique(const PeptideIdentification& identification, PeptideIdentification& filtered_identification);

    /// filter identifications by deviation to the theoretical mass
    static void filterIdentificationsByMzError(const PeptideIdentification& identification, double mass_error, bool unit_ppm, PeptideIdentification& filtered_identification);

    /// only peptides that are in a certain precursor RT range will be kept
    /// Peptides with no RT value will be removed in any case
    static void filterIdentificationsByRT(const std::vector<PeptideIdentification>& identifications, double min_rt, double max_rt, std::vector<PeptideIdentification>& filtered_identifications);

    /// only peptides that are in a certain precursor MZ range will be kept
    /// Peptides with no MZ value will be removed in any case
    static void filterIdentificationsByMZ(const std::vector<PeptideIdentification>& identifications, double min_mz, double max_mz, std::vector<PeptideIdentification>& filtered_identifications);

	  /**
      @brief Filters the peptide hits according to p-values stored in a meta value
    */
    static void filterIdentificationsByPValues(PeptideIdentification& identification, const String& metavalue_key, double p_value = 0.05,
                                               const String& source_tool = "");

	  /**
      @brief Filters the peptide hits according to their predicted RT p-values

      Filters the peptide hits of this ProteinIdentification by the probability (p-value) of a correct ProteinIdentification having a deviation between observed and predicted RT equal to or greater than allowed.
    */
    static void filterIdentificationsByRTPValues(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, double p_value = 0.05);

    /**
      @brief Filters the peptide hits according to their predicted RT p-values of the first dimension (fractionation)

      Filters the peptide hits of this ProteinIdentification by the probability (p-value) of a correct ProteinIdentification having a deviation between observed and predicted RT equal to or greater than allowed.
    */
    static void filterIdentificationsByRTFirstDimPValues(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, double p_value = 0.05);


    /// filters an MS/MS experiment corresponding to the threshold fractions
    template <class PeakT>
    static void filterIdentificationsByThresholds(MSExperiment<PeakT>& experiment, double peptide_threshold_fraction, double protein_threshold_fraction)
    {
      //filter protein hits
      ProteinIdentification temp_protein_identification;
      std::vector<ProteinIdentification> filtered_protein_identifications;

      for (Size j = 0; j < experiment.getProteinIdentifications().size(); j++)
      {
        filterIdentificationsByThreshold(experiment.getProteinIdentifications()[j], protein_threshold_fraction, temp_protein_identification);
        if (!temp_protein_identification.getHits().empty())
        {
          filtered_protein_identifications.push_back(temp_protein_identification);
        }
      }
      experiment.setProteinIdentifications(filtered_protein_identifications);

      //filter peptide hits
      PeptideIdentification temp_identification;
      std::vector<PeptideIdentification> filtered_identifications;

      for (Size i = 0; i < experiment.size(); i++)
      {
        for (Size j = 0; j < experiment[i].getPeptideIdentifications().size(); j++)
        {
          filterIdentificationsByThreshold(experiment[i].getPeptideIdentifications()[j], peptide_threshold_fraction, temp_identification);
          if (!temp_identification.getHits().empty())
          {
            filtered_identifications.push_back(temp_identification);
          }
        }
        experiment[i].setPeptideIdentifications(filtered_identifications);
        filtered_identifications.clear();
      }
    }

    /// filters an MS/MS experiment corresponding to the threshold scores
    template <class PeakT>
    static void filterIdentificationsByScores(MSExperiment<PeakT>& experiment, double peptide_threshold_score, double protein_threshold_score)
    {
      //filter protein hits
      ProteinIdentification temp_protein_identification;
      std::vector<ProteinIdentification> filtered_protein_identifications;

      for (Size j = 0; j < experiment.getProteinIdentifications().size(); j++)
      {
        filterIdentificationsByScore(experiment.getProteinIdentifications()[j], protein_threshold_score, temp_protein_identification);
        if (!temp_protein_identification.getHits().empty())
        {
          filtered_protein_identifications.push_back(temp_protein_identification);
        }
      }
      experiment.setProteinIdentifications(filtered_protein_identifications);

      //filter peptide hits
      PeptideIdentification temp_identification;
      std::vector<PeptideIdentification> filtered_identifications;

      for (Size i = 0; i < experiment.size(); i++)
      {
        for (Size j = 0; j < experiment[i].getPeptideIdentifications().size(); j++)
        {
          filterIdentificationsByScore(experiment[i].getPeptideIdentifications()[j], peptide_threshold_score, temp_identification);
          if (!temp_identification.getHits().empty())
          {
            filtered_identifications.push_back(temp_identification);
          }
        }
        experiment[i].setPeptideIdentifications(filtered_identifications);
        filtered_identifications.clear();
      }
    }

    /// filters an MS/MS experiment corresponding to the best n hits for every spectrum
    template <class PeakT>
    static void filterIdentificationsByBestNHits(MSExperiment<PeakT>& experiment, Size n)
    {
      // @BUG: I don't think it makes sense to remove protein hits and peptide
      // hits independently of each other, since the peptide hits reference the
      // protein hits! (Hendrik)

      //filter protein hits
      ProteinIdentification temp_protein_identification;
      std::vector<ProteinIdentification> filtered_protein_identifications;

      for (Size j = 0; j < experiment.getProteinIdentifications().size(); j++)
      {
        filterIdentificationsByBestNHits(experiment.getProteinIdentifications()[j], n, temp_protein_identification);
        if (!temp_protein_identification.getHits().empty())
        {
          filtered_protein_identifications.push_back(temp_protein_identification);
        }
      }
      experiment.setProteinIdentifications(filtered_protein_identifications);

      //filter peptide hits
      PeptideIdentification temp_identification;
      std::vector<PeptideIdentification> filtered_identifications;

      for (Size i = 0; i < experiment.size(); i++)
      {
        for (Size j = 0; j < experiment[i].getPeptideIdentifications().size(); j++)
        {
          filterIdentificationsByBestNHits(experiment[i].getPeptideIdentifications()[j], n, temp_identification);
          if (!temp_identification.getHits().empty())
          {
            filtered_identifications.push_back(temp_identification);
          }
        }
        experiment[i].setPeptideIdentifications(filtered_identifications);
        filtered_identifications.clear();
      }
    }

    /// filters an MS/MS experiment corresponding to the given proteins
    template <class PeakT>
    static void filterIdentificationsByProteins(MSExperiment<PeakT>& experiment, const std::vector<FASTAFile::FASTAEntry>& proteins)
    {
      std::vector<PeptideIdentification> temp_identifications;
      std::vector<PeptideIdentification> filtered_identifications;
      PeptideIdentification temp_identification;

      for (Size i = 0; i < experiment.size(); i++)
      {
        if (experiment[i].getMSLevel() == 2)
        {
          temp_identifications = experiment[i].getPeptideIdentifications();
          for (Size j = 0; j < temp_identifications.size(); j++)
          {
            filterIdentificationsByProteins(temp_identifications[j], proteins, temp_identification);
            if (!temp_identification.getHits().empty())
            {
              filtered_identifications.push_back(temp_identification);
            }
          }
          experiment[i].setPeptideIdentifications(filtered_identifications);
          filtered_identifications.clear();
        }
      }
    }

    /**
       @brief Update protein groups after protein hits were filtered

       @param groups Input protein groups
       @param hits Available protein hits (all others are removed from the groups)
       @param filtered_groups Output protein groups

       @return Returns whether the groups are still valid (which is the case if only whole groups, if any, were removed).
    */
    static bool updateProteinGroups(
            const std::vector<ProteinIdentification::ProteinGroup>& groups,
            const std::vector<ProteinHit>& hits,
            std::vector<ProteinIdentification::ProteinGroup>& filtered_groups);

  };

} // namespace OpenMS

#endif // OPENMS_FILTERING_ID_IDFILTER_H
