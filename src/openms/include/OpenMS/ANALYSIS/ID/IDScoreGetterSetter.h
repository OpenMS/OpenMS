// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2019.
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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <boost/unordered_map.hpp>

#include <vector>
#include <unordered_set>

namespace OpenMS
{
  /// Used to collect data from the ID structures with the original score as first and
  /// target decoy annotation as second member of the pair. Target = true.
  /// Target+decoy for peptides = target. Protein groups with at least one target = target.
  struct ScoreToTgtDecLabelPairs // Not a typedef to allow forward declaration.
      : public std::vector<std::pair<double, bool>>
  {
    typedef std::vector<std::pair<double, bool>> Base;
    using Base::Base;
  };

  /**
   * @brief A class for extracting and reinserting IDScores from Peptide/ProteinIdentifications and from ConsensusMaps
   */
  class IDScoreGetterSetter
  {

  private:

    template<typename T>
    struct IsIDType
    {
      static bool const value =
          std::is_same<T, PeptideIdentification>::value || std::is_same<T, ProteinIdentification>::value;
    };

    template<typename T>
    struct IsHitType
    {
      static bool const value = std::is_same<T, PeptideHit>::value || std::is_same<T, ProteinHit>::value;
    };

  public:

    /**
     * \defgroup getScoresFunctions Get scores from ID structures for FDR
     * @brief  Fills the scores_labels vector from an ID data structure
     * @param  scores_labels Pairs of scores and boolean target decoy labels to be filled. target = true.
     *
     * Just use the one you need.
     * @{
     */

    //TODO could be done with set of target accessions, too
    //TODO even better: store nr targets and nr decoys when creating the groups!
    //TODO alternative scoring is possible, too (e.g. ratio of tgts vs decoys),
    // this requires templatization of the scores_labels vector though
    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const std::vector<ProteinIdentification::ProteinGroup> &grps,
        const std::unordered_set<std::string> &decoy_accs);

    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const ProteinIdentification &id)
    {

      scores_labels.reserve(scores_labels.size() + id.getHits().size());
      std::transform(id.getHits().begin(), id.getHits().end(),
                     std::back_inserter(scores_labels),
                     [](const ProteinHit &hit)
                     {
                       checkTDAnnotation_(hit);
                       return std::make_pair<double, bool>(hit.getScore(), getTDLabel_(hit));
                     }
      );

    }

    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const PeptideIdentification &id, bool all_hits, int charge, const String &identifier)
    {
      if (id.getIdentifier() == identifier)
      {
        getScores_(scores_labels, id, all_hits, charge);
      }
    }


    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const PeptideIdentification &id, bool all_hits, const String &identifier)
    {
      if (id.getIdentifier() == identifier)
      {
        getScores_(scores_labels, id, all_hits);
      }
    }

    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const PeptideIdentification &id, int charge, const String &identifier)
    {
      if (id.getIdentifier() == identifier)
      {
        getScores_(scores_labels, id, charge);
      }
    }

    template<typename IDType, typename std::enable_if<IsIDType<IDType>::value>::type * = nullptr>
    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const IDType &id, const String &identifier)
    {
      if (id.getIdentifier() == identifier)
      {
        getScores_(scores_labels, id);
      }
    }

    template<class ...Args>
    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const PeptideIdentification &id,
        bool all_hits,
        Args &&... args)
    {
      if (all_hits)
      {
        for (const PeptideHit &hit : id.getHits())
        {
          getScores_(scores_labels, hit, std::forward<Args>(args)...);
        }
      }
      else
      {
        //TODO for speed I assume that they are sorted and first = best.
        //id.sort();
        const PeptideHit &hit = id.getHits()[0];
        getScores_(scores_labels, hit, std::forward<Args>(args)...);
      }
    }

    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const PeptideHit &hit,
        int charge)
    {
      if (charge == hit.getCharge())
      {
        checkTDAnnotation_(hit);
        scores_labels.emplace_back(hit.getScore(), getTDLabel_(hit));
      }
    }

    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const PeptideIdentification &id,
        int charge)
    {
      for (const PeptideHit &hit : id.getHits())
      {
        getScores_(scores_labels, hit, charge);
      }
    }

    template<typename HitType, typename std::enable_if<IsHitType<HitType>::value>::type * = nullptr>
    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const HitType &hit)
    {
      checkTDAnnotation_(hit);
      scores_labels.emplace_back(hit.getScore(), getTDLabel_(hit));
    }

    template<typename IDType, typename std::enable_if<IsIDType<IDType>::value>::type * = nullptr>
    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const IDType &id)
    {
      for (const typename IDType::HitType &hit : id.getHits())
      {
        getScores_(scores_labels, hit);
      }
    }

    template<class ...Args>
    static void getScores_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const std::vector<PeptideIdentification> &ids,
        Args &&... args)
    {
      for (const auto &id : ids)
      {
        getScores_(scores_labels, id, std::forward<Args>(args)...);
      }
    }
    /** @} */

    /**
     * @brief Helper for getting scores in ConsensusMaps
     * @todo allow FeatureMap?
     */
    // GCC-OPT 4.8 -- the following functions can be replaced by a
    // single one with a variadic template, see #4273 and https://gcc.gnu.org/bugzilla/show_bug.cgi?id=41933
    static void getPeptideScoresFromMap_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const ConsensusMap &cmap, bool include_unassigned_peptides)
    {
      std::function<void(const PeptideIdentification &)> f =
          [&](const PeptideIdentification &id) -> void
          { getScores_(scores_labels, id); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void getPeptideScoresFromMap_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const ConsensusMap &cmap, bool include_unassigned_peptides, int charge)
    {
      std::function<void(const PeptideIdentification &)> f =
          [&](const PeptideIdentification &id) -> void
          { getScores_(scores_labels, id, charge); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void getPeptideScoresFromMap_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const ConsensusMap &cmap, bool include_unassigned_peptides, const String &identifier)
    {
      std::function<void(const PeptideIdentification &)> f =
          [&](const PeptideIdentification &id) -> void
          { getScores_(scores_labels, id, identifier); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void getPeptideScoresFromMap_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const ConsensusMap &cmap, bool include_unassigned_peptides, int charge, const String &identifier)
    {
      std::function<void(const PeptideIdentification &)> f =
          [&](const PeptideIdentification &id) -> void
          { getScores_(scores_labels, id, charge, identifier); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void getPeptideScoresFromMap_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const ConsensusMap &cmap, bool include_unassigned_peptides, bool all_hits)
    {
      std::function<void(const PeptideIdentification &)> f =
          [&](const PeptideIdentification &id) -> void
          { getScores_(scores_labels, id, all_hits); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void getPeptideScoresFromMap_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const ConsensusMap &cmap, bool include_unassigned_peptides, bool all_hits, int charge)
    {
      std::function<void(const PeptideIdentification &)> f =
          [&](const PeptideIdentification &id) -> void
          { getScores_(scores_labels, id, all_hits, charge); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void getPeptideScoresFromMap_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const ConsensusMap &cmap, bool include_unassigned_peptides, bool all_hits, const String &identifier)
    {
      std::function<void(const PeptideIdentification &)> f =
          [&](const PeptideIdentification &id) -> void
          { getScores_(scores_labels, id, all_hits, identifier); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void getPeptideScoresFromMap_(
        ScoreToTgtDecLabelPairs &scores_labels,
        const ConsensusMap &cmap, bool include_unassigned_peptides, bool all_hits, int charge, const String &identifier)
    {
      std::function<void(const PeptideIdentification &)> f =
          [&](const PeptideIdentification &id) -> void
          { getScores_(scores_labels, id, all_hits, charge, identifier); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }

    /**
     * @brief For peptide hits, a hit is considered target also if it maps to both
     * a target and a decoy protein (i.e. "target+decoy") as value in the
     * "target_decoy" metavalue e.g. annotated by PeptideIndexer
     */
    static bool getTDLabel_(const MetaInfoInterface &idOrHit)
    {
      return std::string(idOrHit.getMetaValue("target_decoy"))[0] == 't';
    }

    /**
     * \defgroup setScoresFunctions Sets FDRs/qVals
     * @brief  Sets FDRs/qVals from a scores_to_FDR map in the ID data structures
     * @param  scores_to_FDR Maps original score to calculated FDR or q-Value
     * @param  score_type FDR or q-Value
     * @param  higher_better should usually be false @todo remove?
     *
     * Just use the one you need.
     * @{
     */

    template<typename IDType, class ...Args>
    static void setScores_(const std::map<double, double> &scores_to_FDR,
                    std::vector<IDType> &ids,
                    const std::string &score_type,
                    bool higher_better,
                    Args &... args)
    {
      for (auto &id : ids)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better, &args...);
      }
    }

    template<typename IDType>
    static String setScoreType_(IDType &id, const std::string &score_type,
                         bool higher_better)
    {
      String old_score_type = id.getScoreType() + "_score";
      id.setScoreType(score_type);
      id.setHigherScoreBetter(higher_better);
      return old_score_type;
    }

    template<typename IDType>
    static void setScores_(const std::map<double, double> &scores_to_FDR, IDType &id, const std::string &score_type,
                    bool higher_better, bool keep_decoy)
    {
      String old_score_type = setScoreType_(id, score_type, higher_better);

      if (keep_decoy) //in-place set scores
      {
        setScores_(scores_to_FDR, id, old_score_type);
      }
      else
      {
        setScoresAndRemoveDecoys_(scores_to_FDR, id, old_score_type);
      }
    }

    template<typename IDType>
    static void setScores_(const std::map<double, double> &scores_to_FDR, IDType &id,
                    const String &old_score_type)
    {
      std::vector<typename IDType::HitType> &hits = id.getHits();
      for (auto &hit : hits)
      {
        setScore_(scores_to_FDR, hit, old_score_type);
      }
    }

    template<typename IDType, class ...Args>
    static void setScoresAndRemoveDecoys_(const std::map<double, double> &scores_to_FDR, IDType &id,
                                   const String &old_score_type, Args ... args)
    {
      std::vector<typename IDType::HitType> &hits = id.getHits();
      std::vector<typename IDType::HitType> new_hits;
      new_hits.reserve(hits.size());
      for (auto &hit : hits)
      {
        setScoreAndMoveIfTarget_(scores_to_FDR, hit, old_score_type, new_hits, args...);
      }
      hits.swap(new_hits);
    }

    template<typename HitType>
    static void setScore_(const std::map<double, double> &scores_to_FDR, HitType &hit, const std::string &old_score_type)
    {
      hit.setMetaValue(old_score_type, hit.getScore());
      hit.setScore(scores_to_FDR.lower_bound(hit.getScore())->second);
    }

    template<typename IDType>
    static void setScores_(const std::map<double, double> &scores_to_FDR, IDType &id, const std::string &score_type,
                    bool higher_better)
    {
      String old_score_type = setScoreType_(id, score_type, higher_better);
      setScores_(scores_to_FDR, id, old_score_type);
    }

    static void setScores_(const std::map<double, double> &scores_to_FDR,
                    PeptideIdentification &id,
                    const std::string &score_type,
                    bool higher_better,
                    bool keep_decoy,
                    int charge)
    {
      String old_score_type = setScoreType_(id, score_type, higher_better);
      if (keep_decoy) //in-place set scores
      {
        setScores_(scores_to_FDR, id, old_score_type, charge);
      }
      else
      {
        setScoresAndRemoveDecoys_<PeptideIdentification>(scores_to_FDR, id, old_score_type, charge);
      }
    }

    static void setScores_(const std::map<double, double> &scores_to_FDR,
                    PeptideIdentification &id,
                    const std::string &score_type,
                    bool higher_better,
                    bool keep_decoy,
                    int charge,
                    const String &identifier)
    {
      if (id.getIdentifier() == identifier)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better, keep_decoy, charge);
      }
    }

    template<typename IDType>
    static void setScores_(const std::map<double, double> &scores_to_FDR, IDType &id, const std::string &score_type,
                    bool higher_better, bool keep_decoy, const String &identifier)
    {
      if (id.getIdentifier() == identifier)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better, keep_decoy);
      }
    }

    static void setScores_(const std::map<double, double> &scores_to_FDR,
                    PeptideIdentification &id,
                    const std::string &score_type,
                    bool higher_better,
                    int charge,
                    const String &identifier)
    {
      if (id.getIdentifier() == identifier)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better, charge);
      }
    }

    template<typename IDType>
    static void setScores_(const std::map<double, double> &scores_to_FDR, IDType &id, const std::string &score_type,
                    bool higher_better, const String &identifier)
    {
      if (id.getIdentifier() == identifier)
      {
        setScores_(scores_to_FDR, id, score_type, higher_better);
      }
    }

    //TODO could also get a keep_decoy flag when we define what a "decoy group" is -> keep all always for now
    static void setScores_(
        const std::map<double, double> &scores_to_FDR,
        std::vector<ProteinIdentification::ProteinGroup> &grps,
        const std::string &score_type,
        bool higher_better);

    /** @} */

    /**
     * @brief Used when keep_decoy_peptides or proteins is false
     * @tparam HitType ProteinHit or PeptideHit
     * @param scores_to_FDR map from original score to FDR/qVal
     * @param hit The hit (moved to @p new_hits if its a target hit)
     * @param old_score_type to save it in metavalue
     * @param new_hits where to move if target (i.e. target or target+decoy)
     */
    template<typename HitType>
    static void setScoreAndMoveIfTarget_(const std::map<double, double> &scores_to_FDR,
                                  HitType &hit,
                                  const std::string &old_score_type,
                                  std::vector<HitType> &new_hits)
    {
      const String &target_decoy(hit.getMetaValue("target_decoy"));
      if (target_decoy[0] == 't')
      {
        hit.setMetaValue(old_score_type, hit.getScore());
        hit.setScore(scores_to_FDR.lower_bound(hit.getScore())->second);
        new_hits.push_back(std::move(hit));
      } // else do not move over
    }

    /**
    * @brief Used when keep_decoy_peptides is false and charge states are considered
    * @param scores_to_FDR map from original score to FDR/qVal
    * @param hit the PeptideHit itself
    * @param old_score_type to save it in metavalue
    * @param new_hits where to move if target (i.e. target or target+decoy)
    * @param charge If only peptides with charge X are currently considered
    */
    static void setScoreAndMoveIfTarget_(const std::map<double, double> &scores_to_FDR,
                                  PeptideHit &hit,
                                  const std::string &old_score_type,
                                  std::vector<PeptideHit> &new_hits,
                                  int charge)
    {
      if (charge == hit.getCharge())
      {
        const String &target_decoy(hit.getMetaValue("target_decoy"));
        if (target_decoy[0] == 't')
        {
          hit.setMetaValue(old_score_type, hit.getScore());
          hit.setScore(scores_to_FDR.lower_bound(hit.getScore())->second);
          new_hits.push_back(std::move(hit));
        } // else do not move over
      }
      else // different charge, move over unchanged to process later at correct charge.
      {
        new_hits.push_back(std::move(hit));
      }
    }

    /**
     * @brief Helper for applying set Scores on ConsensusMaps
     * @tparam Args optional additional arguments (charge, run ID)
     * @param scores_to_FDR maps original scores to FDR
     * @param cmap the ConsensusMap
     * @param score_type FDR or q-Value
     * @param higher_better usually false
     * @param keep_decoy read from Param object
     * @param args optional additional arguments (int charge, string run ID)
     */
    // GCC-OPT 4.8 -- the following functions can be replaced by a
    // single one with a variadic template, see #4273 and https://gcc.gnu.org/bugzilla/show_bug.cgi?id=41933
    static void setPeptideScoresForMap_(const std::map<double, double> &scores_to_FDR,
                                 ConsensusMap &cmap,
                                 bool include_unassigned_peptides,
                                 const std::string &score_type,
                                 bool higher_better,
                                 bool keep_decoy)
    {
      //Note: Gcc4.8 cannot handle variadic templates in lambdas
      std::function<void(PeptideIdentification &)> f =
          [&](PeptideIdentification &id) -> void
          { setScores_(scores_to_FDR, id, score_type,
              higher_better, keep_decoy); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void setPeptideScoresForMap_(const std::map<double, double> &scores_to_FDR,
                                        ConsensusMap &cmap,
                                        bool include_unassigned_peptides,
                                        const std::string &score_type,
                                        bool higher_better,
                                        bool keep_decoy,
                                        int charge)
    {
      //Note: Gcc4.8 cannot handle variadic templates in lambdas
      std::function<void(PeptideIdentification &)> f =
          [&](PeptideIdentification &id) -> void
          { setScores_(scores_to_FDR, id, score_type,
                       higher_better, keep_decoy, charge); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void setPeptideScoresForMap_(const std::map<double, double> &scores_to_FDR,
                                        ConsensusMap &cmap,
                                        bool include_unassigned_peptides,
                                        const std::string &score_type,
                                        bool higher_better,
                                        bool keep_decoy,
                                        const String& run_identifier)
    {
      //Note: Gcc4.8 cannot handle variadic templates in lambdas
      std::function<void(PeptideIdentification &)> f =
          [&](PeptideIdentification &id) -> void
          { setScores_(scores_to_FDR, id, score_type,
                       higher_better, keep_decoy, run_identifier); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }
    static void setPeptideScoresForMap_(const std::map<double, double> &scores_to_FDR,
                                        ConsensusMap &cmap,
                                        bool include_unassigned_peptides,
                                        const std::string &score_type,
                                        bool higher_better,
                                        bool keep_decoy,
                                        int charge,
                                        const String& run_identifier)
    {
      //Note: Gcc4.8 cannot handle variadic templates in lambdas
      std::function<void(PeptideIdentification &)> f =
          [&](PeptideIdentification &id) -> void
          { setScores_(scores_to_FDR, id, score_type,
                       higher_better, keep_decoy, charge, run_identifier); };
      cmap.applyFunctionOnPeptideIDs(f, include_unassigned_peptides);
    }

    /**
     * @brief To check the metavalues before we do anything
     * @param id_or_hit Any Object with MetaInfoInterface. Specifically ID or Hit Type here.
     * @throws Exception::MissingInformation if target_decoy annotation does not exist
     */
    static void checkTDAnnotation_(const MetaInfoInterface &id_or_hit)
    {
      if (!id_or_hit.metaValueExists("target_decoy"))
      {
        throw Exception::MissingInformation(__FILE__,
                                            __LINE__,
                                            OPENMS_PRETTY_FUNCTION,
                                            "Meta value 'target_decoy' does not exist in all ProteinHits! Reindex the idXML file with 'PeptideIndexer'");
      }
    }
  };
} // namespace OpenMS
