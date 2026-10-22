// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/IdentificationDataEdits.h>

#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
  using ID = IdentificationData;

  namespace
  {
    /// The run of @p data with @p uuid
    ID::Run& runOf(ID& data, const std::string& uuid)
    {
      auto* run = data.findRunByUuid(uuid);
      if (run == nullptr)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The run is no longer part of the identification data", uuid);
      }
      return *run;
    }
  } // namespace

  IdentificationDataEdits::IdentificationDataEdits(IdentificationData& data) : data_(data), results_(data.getInferenceResults())
  {
  }

  ID::RunSettings& IdentificationDataEdits::settings(const ID::Run& run)
  {
    return settings_.try_emplace(run.getUuid(), run.getSettings()).first->second;
  }

  std::vector<ID::DatabaseSequence>& IdentificationDataEdits::databaseSequences(const ID::Run& run)
  {
    if (!run.getDatabaseSequences())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The run has no database sequences", run.getIdentifier());
    }
    return sequences_.try_emplace(run.getUuid(), *run.getDatabaseSequences()).first->second;
  }

  ID::Observation& IdentificationDataEdits::observation(const ID::Run& run, const ID::Identification& query)
  {
    return observations_.try_emplace({run.getUuid(), query.getId().value}, query.getObservation()).first->second;
  }

  ID::MatchData& IdentificationDataEdits::match(const ID::Run& run, const ID::Match& match)
  {
    return matches_.try_emplace({run.getUuid(), match.getId().value}, match.getData()).first->second;
  }

  bool IdentificationDataEdits::commit()
  {
    // the edits go to a copy, which replaces the data when all are valid
    ID edited(data_);
    bool changed = false;
    for (const auto& [uuid, settings] : settings_)
    {
      auto& run = runOf(edited, uuid);
      if (run.getSettings() == settings) continue;
      run.setSettings(settings);
      changed = true;
    }
    for (const auto& [uuid, sequences] : sequences_)
    {
      auto& run = runOf(edited, uuid);
      if (run.getDatabaseSequences() && *run.getDatabaseSequences() == sequences) continue;
      run.setDatabaseSequences(sequences);
      changed = true;
    }
    for (const auto& [key, observation] : observations_)
    {
      auto& run = runOf(edited, key.first);
      const ID::QueryId query {key.second};
      if (run.getIdentification(query).getObservation() == observation) continue;
      run.replaceObservation(query, observation);
      changed = true;
    }
    for (const auto& [key, data] : matches_)
    {
      auto& run = runOf(edited, key.first);
      const ID::MatchId match {key.second};
      if (run.getMatch(match).getData() == data) continue;
      run.replaceMatch(match, data);
      changed = true;
    }
    // a protein score definition no longer applies to proteins of another score type (see IdentificationDataAdapter::editInferenceResults())
    auto results = results_;
    for (auto& result : results)
    {
      if (result.protein_score && result.protein_score->name != result.proteins.getScoreType()) result.protein_score.reset();
    }
    if (edited.getInferenceResults() != results)
    {
      edited.clearInferenceResults();
      for (auto& result : results) edited.addInferenceResult(std::move(result));
      changed = true;
    }
    if (!changed) return false;
    edited.validate();
    data_.swap(edited);
    return true;
  }
} // namespace OpenMS
