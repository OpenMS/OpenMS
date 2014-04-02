// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Authors: Nico Pfeifer, Mathias Walzer$
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_ID_IDFILTER_H
#define OPENMS_FILTERING_ID_IDFILTER_H

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <vector>
#include <climits>

namespace OpenMS
{
  /**
    @brief used to filter identifications by different criteria

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

    /// gets the best scoring peptide hit from a vector of peptide identifications
    /// @param identifications Vector of peptide ids, each containing one or more peptide hits
    /// @param assume_sorted are hits sorted by score (best score first) already? This allows for faster query, since only the first hit needs to be looked at
    /// @return true if a hit was present, false otherwise
    template <class IdentificationType>
    bool getBestHit(const std::vector<IdentificationType> identifications, bool assume_sorted, PeptideHit& best_hit)
    {
      if (identifications.size() == 0) return false;

      bool is_higher_score_better = identifications[0].isHigherScoreBetter();
      double best_score = (is_higher_score_better ? -1 : 1) * std::numeric_limits<double>::max(); // worst score we can think of

      Size best_i_index(0), best_h_index(0);
      Size max_h(-1);
      // determine best scoring hit
      for (Size i = 0; i != identifications.size(); ++i)
      {
        if (identifications[i].getHits().size() == 0) continue; // empty hits
        
        is_higher_score_better = identifications[i].isHigherScoreBetter();
        max_h = (assume_sorted ? 1 : identifications[i].getHits().size());        
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

      if (max_h == -1) return false;// all hits were empty 

      best_hit = identifications[best_i_index].getHits()[best_h_index];
      return true;

    }

    /// filters a ProteinIdentification or PeptideIdentification by only allowing peptides/proteins which reach a score above @p threshold_fraction * SignificanceThreshold
    template <class IdentificationType>
    void filterIdentificationsByThreshold(const IdentificationType& identification, double threshold_fraction, IdentificationType& filtered_identification)
    {
      typedef typename IdentificationType::HitType HitType;
      std::vector<HitType> temp_hits;
      std::vector<HitType> filtered_hits;

      filtered_identification = identification;
      filtered_identification.setHits(std::vector<HitType>());

      for (typename std::vector<HitType>::const_iterator it = identification.getHits().begin();
           it != identification.getHits().end();
           ++it)
      {
        if (it->getScore() >= threshold_fraction * identification.getSignificanceThreshold())
        {
          filtered_hits.push_back(*it);
        }
      }

      if (!filtered_hits.empty())
      {
        filtered_identification.setHits(filtered_hits);
        filtered_identification.assignRanks();
      }
    }

    /**
      @brief filters a ProteinIdentification or PeptideIdentification corresponding to the @p threshold_score

      If the method higherScoreBetter() returns true for the IdentificationType all hits with a score
      smaller than @p threshold_score are removed. Otherwise all hits with a score bigger than
      @p threshold_score are removed.
    */
    template <class IdentificationType>
    void filterIdentificationsByScore(const IdentificationType& identification, double threshold_score, IdentificationType& filtered_identification)
    {
      typedef typename IdentificationType::HitType HitType;
      std::vector<HitType> temp_hits;
      std::vector<HitType> filtered_hits;

      filtered_identification = identification;
      filtered_identification.setHits(std::vector<HitType>());

      for (typename std::vector<HitType>::const_iterator it = identification.getHits().begin();
           it != identification.getHits().end();
           ++it)
      {
        if (identification.isHigherScoreBetter())
        {
          if (it->getScore() >= threshold_score)
          {
            filtered_hits.push_back(*it);
          }
        }
        else
        {
          if (it->getScore() <= threshold_score)
          {
            filtered_hits.push_back(*it);
          }
        }
      }

      if (!filtered_hits.empty())
      {
        filtered_identification.setHits(filtered_hits);
        filtered_identification.assignRanks();
      }
    }

    /**
      @brief filters a ProteinIdentification or PeptideIdentification corresponding to the score.

      If the method higherScoreBetter() returns true for the IdentificationType the
      n highestscoring hits are kept. Otherwise the n lowest scoring hits are kept.
    */
    template <class IdentificationType>
    void filterIdentificationsByBestNHits(const IdentificationType& identification, Size n, IdentificationType& filtered_identification)
    {
      typedef typename IdentificationType::HitType HitType;
      std::vector<HitType> temp_hits;
      std::vector<HitType> filtered_hits;
      Size count = 0;

      IdentificationType temp_identification = identification;
      temp_identification.sort(); // .. by score

      filtered_identification = identification;
      filtered_identification.setHits(std::vector<HitType>());


      typename std::vector<HitType>::const_iterator it = temp_identification.getHits().begin();
      while (it != temp_identification.getHits().end()
            && count < n)
      {
        filtered_hits.push_back(*it);
        ++it;
        ++count;
      }

      if (!filtered_hits.empty())
      {
        filtered_identification.setHits(filtered_hits);
        filtered_identification.assignRanks();
      }
    }

    /**
      @brief filters a ProteinIdentification or PeptideIdentification corresponding to the score.

      If the method higherScoreBetter() returns true for the IdentificationType the
      n to m highestscoring hits are kept. Otherwise the n to m lowest scoring hits are kept.
      This method is useful if a range of higher hits are used for decoy fairness analysis.
    */
    template <class IdentificationType>
    void filterIdentificationsByBestNToMHits(const IdentificationType& identification, Size n, Size m, IdentificationType& filtered_identification)
    {
      if (n > m)
      {
        std::swap(n, m);
      }

      typedef typename IdentificationType::HitType HitType;
      std::vector<HitType> filtered_hits;

      IdentificationType temp_identification = identification;
      temp_identification.sort(); // .. by score

      filtered_identification = identification;
      filtered_identification.setHits(std::vector<HitType>());

      const std::vector<HitType>& hits = temp_identification.getHits();
      for (Size i = n - 1; n <= m - 1; ++i)
      {
        if (i >= hits.size())
        {
          break;
        }
        filtered_hits.push_back(hits[i]);
      }

      if (!filtered_hits.empty())
      {
        filtered_identification.setHits(filtered_hits);
        filtered_identification.assignRanks();
      }
    }

    /// filters a PeptideIdentification keeping only the best scoring hits (if strict is set, keeping only the best hit only if it is the only hit with that score)
    void filterIdentificationsByBestHits(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, bool strict = false);

    /**
       @brief Checks whether a meta value of the peptide identification is within a given range

       Useful for filtering by precursor RT or m/z.

       @param identification The peptide ID to check
       @param key Key (name) for the meta value
       @param low Lower boundary (inclusive)
       @param high Upper boundary (inclusive)
       @param missing What to return when the meta value is missing

       @returns Whether the peptide ID passes the check
    */
    bool filterIdentificationsByMetaValueRange(const PeptideIdentification& identification, const String& key, double low, double high, bool missing = false);
    
    /// filters a PeptideIdentification corresponding to the given proteins
    /// PeptideHits with no matching @em proteins are removed.
    /// Matching is done either based on accessions or on sequence (if no accessions are given, or @em no_protein_identifiers is set)
    void filterIdentificationsByProteins(const PeptideIdentification& identification, const std::vector<FASTAFile::FASTAEntry>& proteins, PeptideIdentification& filtered_identification, bool no_protein_identifiers = false);

    /// filters a ProteinIdentification corresponding to the given proteins
    /// ProteinHits with no matching @em proteins are removed.
    /// Matching is done based on accessions only
    void filterIdentificationsByProteins(const ProteinIdentification& identification, const std::vector<FASTAFile::FASTAEntry>& proteins, ProteinIdentification& filtered_identification);

    /// removes all peptide hits having a sequence equal to a element in <code>peptides</code>
    void filterIdentificationsByExclusionPeptides(const PeptideIdentification& identification, const std::set<String>& peptides, PeptideIdentification& filtered_identification);

    /// only peptides having a length equal to or greater than 'min_length' will be kept, if 'max_length' is actually greater or equal 'min_length', also only peptides having a length less or equal to 'max_length' will be kept
    void filterIdentificationsByLength(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, Size min_length, Size max_length = UINT_MAX);

    /// only peptides that have a charge equal to or greater than 'charge' will be kept
    void filterIdentificationsByCharge(const PeptideIdentification& identification, Int charge, PeptideIdentification& filtered_identification);

    /// only peptides having a variable modification will be kept
    void filterIdentificationsByVariableModifications(const PeptideIdentification& identification, const std::vector<String>& fixed_modifications, PeptideIdentification& filtered_identification);

    /// only protein hits in 'identification' which are referenced by a peptide in 'peptide_identifications' are kept
    void removeUnreferencedProteinHits(const ProteinIdentification& identification, const std::vector<PeptideIdentification> peptide_identifications, ProteinIdentification& filtered_identification);

    /// only peptide hits in 'peptide_identifications' which are referenced by a protein in 'identification' are kept
    void removeUnreferencedPeptideHits(const ProteinIdentification& identification, std::vector<PeptideIdentification>& peptide_identifications, bool delete_unreferenced_peptide_hits = false);

    /// if a peptide hit occurs more than once per PSM, only one instance is kept
    void filterIdentificationsUnique(const PeptideIdentification& identification, PeptideIdentification& filtered_identification);

    /// filter identifications by deviation to the theoretical mass
    void filterIdentificationsByMzError(const PeptideIdentification& identification, double mass_error, bool unit_ppm, PeptideIdentification& filtered_identification);

    /**
          @brief Filters the peptide hits according to their predicted rt p-values

          Filters the peptide hits of this ProteinIdentification by the
          probability (p-value) of a correct ProteinIdentification having a deviation between
          observed and predicted rt equal or bigger than allowed.
      */
    void filterIdentificationsByRTPValues(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, double p_value = 0.05);

    /**
          @brief Filters the peptide hits according to their predicted rt p-values of the first dimension

          Filters the peptide hits of this ProteinIdentification by the
          probability (p-value) of a correct ProteinIdentification having a deviation between
          observed and predicted rt equal or bigger than allowed.
      */
    void filterIdentificationsByRTFirstDimPValues(const PeptideIdentification& identification, PeptideIdentification& filtered_identification, double p_value = 0.05);

    /// filters an MS/MS experiment corresponding to the threshold_fractions
    template <class PeakT>
    void filterIdentificationsByThresholds(MSExperiment<PeakT>& experiment, double peptide_threshold_fraction, double protein_threshold_fraction)
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

    /// filters an MS/MS experiment corresponding to the threshold_fractions
    template <class PeakT>
    void filterIdentificationsByScores(MSExperiment<PeakT>& experiment, double peptide_threshold_score, double protein_threshold_score)
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
    void filterIdentificationsByBestNHits(MSExperiment<PeakT>& experiment, Size n)
    {
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
    void filterIdentificationsByProteins(MSExperiment<PeakT>& experiment, const std::vector<FASTAFile::FASTAEntry>& proteins)
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

  };

} // namespace OpenMS

#endif // OPENMS_FILTERING_ID_IDFILTER_H
