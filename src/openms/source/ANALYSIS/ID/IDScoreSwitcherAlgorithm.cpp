// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDScoreSwitcherAlgorithm.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ID/IdentificationDataAdapter.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <unordered_map>

using namespace std;
namespace OpenMS
{

  IDScoreSwitcherAlgorithm::IDScoreSwitcherAlgorithm() :
      IDScoreSwitcherAlgorithm::DefaultParamHandler("IDScoreSwitcherAlgorithm")
  {
    defaults_.setValue("new_score", "", "Name of the meta value to use as the new score");
    defaults_.setValue("new_score_orientation", "", "Orientation of the new score (are higher or lower values better?)");
    defaults_.setValidStrings("new_score_orientation", {"lower_better","higher_better"});
    defaults_.setValue("new_score_type", "", "Name to use as the type of the new score (default: same as 'new_score')");
    defaults_.setValue("old_score", "", "Name to use for the meta value storing the old score (default: old score type)");
    defaults_.setValue("proteins", "false", "Apply to protein scores instead of PSM scores");
    defaults_.setValidStrings("proteins", {"true","false"});
    defaultsToParam_();
    updateMembers_();
  }

  void IDScoreSwitcherAlgorithm::updateMembers_()
  {
    new_score_ = param_.getValue("new_score").toString();
    new_score_type_ = param_.getValue("new_score_type").toString();
    old_score_ = param_.getValue("old_score").toString();
    higher_better_ = (param_.getValue("new_score_orientation").toString() ==
                      "higher_better");

    if (new_score_type_.empty()) new_score_type_ = new_score_;
  }

  namespace
  {
    using ID = IdentificationData;

    /// Run @p edit on the identification data of @p cmap (see IdentificationDataConverter::editAsIdentificationData())
    template<class Edit>
    void editIdentificationData(ConsensusMap& cmap, Edit&& edit)
    {
      IdentificationDataConverter::editAsIdentificationData(cmap, [&](ConsensusMap& map) { edit(map.getIdentificationData()); });
    }

    /// The first identification that a feature of @p cmap links, with its first linked match
    std::optional<ID::QueryMatches> firstAssigned(const ConsensusMap& cmap)
    {
      for (const auto& feature : cmap)
      {
        auto linked = feature.getLinkedIdentifications(cmap.getIdentificationData());
        if (! linked.empty()) return std::move(linked.front());
      }
      return std::nullopt;
    }
  } // namespace

  void IDScoreSwitcherAlgorithm::switchToGeneralScoreType(ConsensusMap& cmap, ScoreType type, Size& counter, bool /* unassigned_peptides_too */)
  {
    editIdentificationData(cmap, [&](ID& data) {
      const auto primary = data.getPrimaryScoreDefinition();
      std::string new_type;
      // The first identification of a feature, as a peptide identification (with its scores as meta values), tells
      // which score to switch to; features are tried in order until one has it.
      for (const auto& feature : cmap)
      {
        if (! primary) break;
        const auto linked = feature.getLinkedIdentifications(data);
        if (linked.empty() || linked.front().matches.empty()) continue;
        const auto& run = *linked.front().run;
        const auto& match = *linked.front().matches.front();
        PeptideIdentification id;
        id.setScoreType(primary->name);
        id.setHigherScoreBetter(primary->higher_better);
        auto hit = IdentificationDataAdapter::materializePeptide(run, match, *run.getPrimaryScore());
        const auto scores = run.getScores(match);
        for (Size i = 0; i < scores.size(); ++i)
          if (scores[i] && ! (run.getScoreDefinitions()[i] == *primary)) hit.setMetaValue(run.getScoreDefinitions()[i].name, *scores[i]);
        id.insertHit(hit);
        const auto sr = findScoreType(id, type);
        if (sr.is_main_score_type) return;
        if (! sr.score_name.empty())
        {
          new_type = sr.score_name;
          break;
        }
      }
      if (new_type.empty())
      {
        std::string msg = "First encountered ID does not have the requested score type.";
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg);
      }
      new_score_type_ = StringUtils::hasSuffix(new_type, "_score") ? StringUtils::chop(new_type, 6) : new_type;
      new_score_ = new_type;
      if (higher_better_ != Scores::isHigherBetter(type))
      {
        OPENMS_LOG_WARN << "Requested score type does not match the expected score direction. Correcting!\n";
        higher_better_ = Scores::isHigherBetter(type);
      }
      switchScores_(data, counter);
    });
  }

  void IDScoreSwitcherAlgorithm::switchScores(ConsensusMap& cmap, Size& counter, bool /* unassigned_peptides_too */)
  {
    editIdentificationData(cmap, [&](ID& data) {
      const auto primary = data.getPrimaryScoreDefinition();
      if (! primary || primary->name == new_score_) return; // correct score or category already set
      switchScores_(data, counter);
    });
  }

  void IDScoreSwitcherAlgorithm::switchScores_(IdentificationData& data, Size& counter) const
  {
    // The new score is the meta value new_score_ of every match (or the score of that name); the previous one becomes
    // the meta value old_score_ (default: its name), as switchScores() of peptide identifications does.
    ID::ScoreDefinition definition;
    definition.name = new_score_type_;
    definition.higher_better = higher_better_;
    IdentificationDataAdapter::replacePrimaryScore(data, definition,
      [&](const ID::Run& run, const ID::Match& match, double) {
        ++counter;
        const auto& definitions = run.getScoreDefinitions();
        for (Size i = 0; i < definitions.size(); ++i)
          if (definitions[i].name == new_score_)
          {
            const auto value = run.getScores(match)[i];
            if (value) return *value;
          }
        if (! match.metaValueExists(new_score_))
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Meta value '" + new_score_ + "' not found for " + match.representation);
        }
        return double(match.getMetaValue(new_score_));
      },
      "", true, old_score_);
  }

  void IDScoreSwitcherAlgorithm::determineScoreNameOrientationAndType(const ConsensusMap& cmap, std::string& name, bool& higher_better,
                                                                       ScoreType& score_type, bool include_unassigned)
  {
    name = "";
    higher_better = true;
    std::optional<ConsensusMap> converted;
    const ConsensusMap& map = IdentificationDataConverter::withIdentificationData(cmap, converted);
    const auto primary = map.getIdentificationData().getPrimaryScoreDefinition();
    // The main score of the identification data is that of the assigned identifications, and of the unassigned ones.
    if (! primary || (! firstAssigned(map) && (! include_unassigned || map.getUnassignedIdentifications().empty()))) return;
    name = primary->name;
    higher_better = primary->higher_better;
    // look up the score category ("RAW", "PEP", "q-value", etc.) for the given score name
    Scores::findIDTypeByName(name, score_type);
  }

  std::vector<std::string> IDScoreSwitcherAlgorithm::getScoreNames()
  {
    return Scores::getAllIDScoreNames();
  }


} // namespace OpenMS
