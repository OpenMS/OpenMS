// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledIdentifications.h"

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FEATUREFINDER/MultiplexResolverAlgorithm.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace OpenMS
{
namespace
{
using ID = IdentificationData;

/// The FAIMS CV a multiplet, feature or identification was acquired at (none for conventional data)
std::optional<double> compensationVoltage(const MetaInfoInterface& value)
{
  if (!value.metaValueExists(Constants::UserParam::FAIMS_CV)) { return std::nullopt; }
  return static_cast<double>(value.getMetaValue(Constants::UserParam::FAIMS_CV));
}
} // namespace

void MS1LabeledIdentifications::moveIdentifications(ConsensusMap& multiplets, const std::vector<Move>& moves)
{
  if (moves.empty()) { return; }
  ID& data = multiplets.getIdentificationData();
  struct Copy
  {
    std::string run;
    ID::SourceId source;
    ID::Observation observation;
    std::vector<std::pair<ID::MatchData, std::vector<std::optional<double>>>> matches;
    std::optional<Size> multiplet;
  };
  std::vector<Copy> copies;
  std::set<ID::QueryReference> originals;
  for (const Move& move : moves)
  {
    const ID::Run* run = data.findRunByUuid(move.identification.run_uuid);
    const ID::Identification* query = run == nullptr ? nullptr : run->findIdentification(move.identification.query);
    if (query == nullptr || !originals.insert(move.identification).second)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "An identification to move is unknown or moved twice.",
                                    move.identification.run_uuid);
    }
    if (move.multiplet && *move.multiplet >= multiplets.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "An identification is to be moved to a multiplet that does not exist.",
                                    StringUtils::toStr(*move.multiplet));
    }
    for (const auto& source : run->getSources())
    {
      if (std::none_of(source.identifications.begin(), source.identifications.end(),
                       [&](const ID::Identification& candidate) { return candidate.getId() == query->getId(); }))
      {
        continue;
      }
      Copy copy {run->getIdentifier(), source.id, query->getObservation(), {}, move.multiplet};
      if (move.map_index) { copy.observation.setMetaValue("map_index", *move.map_index); }
      else { copy.observation.removeMetaValue("map_index"); }
      for (const auto& match : query->getMatches()) { copy.matches.emplace_back(match.getData(), run->getScores(match)); }
      copies.push_back(std::move(copy));
      break;
    }
  }
  UInt64 next = 1;
  for (const auto& run : data.getRuns()) { next = std::max(next, run.getNextQueryId()); }
  for (auto& copy : copies)
  {
    auto& run = data.getRun(copy.run);
    const ID::QueryId query = run.importIdentification(copy.source, ID::QueryId {next++}, std::move(copy.observation));
    std::vector<ID::MatchId> matches;
    for (const auto& [match, scores] : copy.matches) { matches.push_back(run.addMatch(query, match, scores)); }
    if (!copy.multiplet) { continue; }
    ConsensusFeature& multiplet = multiplets[*copy.multiplet];
    multiplet.addIDQuery({run.getUuid(), query});
    for (const ID::MatchId& match : matches) { multiplet.addIDMatch({run.getUuid(), match}); }
  }
  multiplets.eraseIdentifications([&](const ID::Run& run, const ID::Identification& query) { return originals.contains({run.getUuid(), query.getId()}); });
}

Size MS1LabeledIdentifications::attachToUndetectedChannels(ConsensusMap& multiplets, const FeatureMap& channels, const MultiplexResolverAlgorithm& resolver,
                                                           Size multiplicity, const Tolerances& tolerances)
{
  std::unordered_map<UInt64, Size> multiplet_of;
  for (Size m = 0; m < multiplets.size(); ++m)
  {
    for (const FeatureHandle& handle : multiplets[m].getFeatures()) { multiplet_of.emplace(handle.getUniqueId(), m); }
  }
  std::vector<std::pair<double, const Feature*>> by_mz;
  by_mz.reserve(channels.size());
  for (const Feature& channel : channels) { by_mz.emplace_back(channel.getMZ(), &channel); }
  std::sort(by_mz.begin(), by_mz.end(), [](const auto& a, const auto& b) { return a.first < b.first; });

  auto tolerance_at = [&](double mz) { return tolerances.ppm ? mz * tolerances.mz * 1e-6 : tolerances.mz; };

  std::vector<Move> moves;
  for (const auto& entry : multiplets.getUnassignedIdentifications())
  {
    if (entry.matches.empty() || entry.matches[0]->charge == 0 || !entry.query->rt) { continue; }
    const ID::Match& match = *entry.matches[0];
    const double rt = *entry.query->rt;
    Size own = 0;
    const std::vector<double> mzs = resolver.expectedChannelMZs(AASequence::fromString(match.representation), match.charge, own);
    const std::optional<double> cv = compensationVoltage(*entry.query);

    // the closest channel feature at the m/z of the peptide in another channel
    const Feature* best = nullptr;
    double best_error = std::numeric_limits<double>::max();
    for (Size j = 0; j < mzs.size(); ++j)
    {
      if (j == own) { continue; }
      const double tolerance = tolerance_at(mzs[j]);
      auto it = std::lower_bound(by_mz.begin(), by_mz.end(), mzs[j] - tolerance, [](const auto& a, double mz) { return a.first < mz; });
      for (; it != by_mz.end() && it->first <= mzs[j] + tolerance; ++it)
      {
        const Feature& channel = *it->second;
        if (channel.getCharge() != match.charge || compensationVoltage(channel) != cv) { continue; }
        const DBoundingBox<2> box = channel.getConvexHull().getBoundingBox();
        if (rt < box.minX() - tolerances.rt || rt > box.maxX() + tolerances.rt) { continue; }
        // in a complete multiplet the identification's channel was detected and did not match it (placeholders added
        // below for earlier identifications have no unique id and do not count)
        const auto& handles = multiplets[multiplet_of.at(channel.getUniqueId())].getFeatures();
        if (static_cast<Size>(std::count_if(handles.begin(), handles.end(), [](const FeatureHandle& h) { return h.hasValidUniqueId(); })) >= multiplicity)
        {
          continue;
        }
        const double error = std::abs(channel.getMZ() - mzs[j]);
        if (error < best_error)
        {
          best = &channel;
          best_error = error;
        }
      }
    }
    if (best == nullptr) { continue; }

    // the identification's own channel joins the multiplet (unless an earlier identification added it already)
    const Size m = multiplet_of.at(best->getUniqueId());
    ConsensusFeature& multiplet = multiplets[m];
    const double own_mz = mzs[own];
    std::vector<FeatureHandle> handles(multiplet.getFeatures().begin(), multiplet.getFeatures().end());
    if (std::none_of(handles.begin(), handles.end(), [&](const FeatureHandle& h) { return std::abs(h.getMZ() - own_mz) <= tolerance_at(own_mz); }))
    {
      FeatureHandle placeholder;
      placeholder.setMZ(own_mz);
      placeholder.setRT(multiplet.getRT());
      placeholder.setCharge(match.charge);
      placeholder.setIntensity(0.0);
      handles.push_back(placeholder);
      // channels in m/z order, numbered from 0 like a detected (partial) pattern; the lightest is the anchor
      std::sort(handles.begin(), handles.end(), [](const FeatureHandle& a, const FeatureHandle& b) { return a.getMZ() < b.getMZ(); });
      ConsensusFeature::HandleSetType renumbered;
      for (Size i = 0; i < handles.size(); ++i)
      {
        handles[i].setMapIndex(i);
        renumbered.insert(handles[i]);
      }
      multiplet.setFeatures(std::move(renumbered));
      multiplet.setMZ(handles.front().getMZ());
    }
    // the channel is settled by keepMultiAssignedOnClosest()
    moves.push_back({{entry.run->getUuid(), entry.query->getId()}, m, 0});
  }
  // the identifications join their multiplets after the ones mapped onto them
  moveIdentifications(multiplets, moves);
  return moves.size();
}

Size MS1LabeledIdentifications::keepMultiAssignedOnClosest(ConsensusMap& multiplets)
{
  struct Candidate
  {
    ID::QueryReference identification;
    double mz_ppm;
    double rt_diff;
  };
  std::map<std::string, std::vector<Candidate>> by_spectrum_match; // (spectrum, peptidoform, charge) -> where it sits

  // Every assignment is an identification of its own (IDMapper with 'annotate_ids_with_subelements').
  ID& data = multiplets.getIdentificationData();
  for (const ConsensusFeature& cf : multiplets)
  {
    for (const auto& entry : cf.getLinkedIdentifications(data))
    {
      if (entry.matches.empty()) { continue; }
      const ID::Identification& id = *entry.query;
      const double id_rt = id.rt.value_or(0.0);

      // the channel closest to the identified peptidoform is the one the spectrum belongs to
      const ID::Match& first = *entry.matches[0];
      const double id_mz = first.charge != 0 ? AASequence::fromString(first.representation).getMZ(first.charge) : id.mz.value_or(0.0);
      double best_ppm = std::numeric_limits<double>::infinity(), best_rt = 0.0;
      Size best_map_index = 0;
      for (const FeatureHandle& handle : cf.getFeatures())
      {
        const double ppm = std::abs(handle.getMZ() - id_mz) / handle.getMZ() * 1e6;
        if (ppm < best_ppm)
        {
          best_ppm = ppm;
          best_rt = std::abs(handle.getRT() - id_rt);
          best_map_index = handle.getMapIndex();
        }
      }
      if (!id.data_id.empty())
      {
        const std::string key = id.data_id + '\t' + first.representation + '\t' + StringUtils::toStr(first.charge);
        by_spectrum_match[key].push_back({{entry.run->getUuid(), id.getId()}, best_ppm, best_rt});
      }
      if (!cf.getFeatures().empty())
      {
        ID::Observation observation = id;
        observation.setMetaValue("map_index", best_map_index);
        data.getRun(entry.run->getIdentifier()).replaceObservation(id.getId(), observation);
      }
    }
  }

  std::set<ID::QueryReference> to_remove; // the copies to drop
  Size shared = 0;
  for (auto& [key, candidates] : by_spectrum_match)
  {
    if (candidates.size() < 2) { continue; }
    ++shared;
    const auto closest = std::min_element(candidates.begin(), candidates.end(),
                                          [](const Candidate& a, const Candidate& b)
                                          { return std::tie(a.mz_ppm, a.rt_diff) < std::tie(b.mz_ppm, b.rt_diff); });
    for (auto it = candidates.begin(); it != candidates.end(); ++it)
    {
      if (it != closest) { to_remove.insert(it->identification); }
    }
  }
  if (!to_remove.empty())
  {
    multiplets.eraseIdentifications([&](const ID::Run& run, const ID::Identification& id) { return to_remove.contains({run.getUuid(), id.getId()}); });
  }
  return shared;
}

Size MS1LabeledIdentifications::keepConsistentIdentifications(ConsensusMap& multiplets, const MultiplexResolverAlgorithm& resolver)
{
  const ID& data = multiplets.getIdentificationData();
  std::vector<Move> moves;
  for (const auto& cf : multiplets)
  {
    const auto ids = cf.getLinkedIdentifications(data);
    if (ids.size() < 2) { continue; }
    const std::vector<bool> consistent = resolver.checkIdentifications(cf, data);
    if (std::none_of(consistent.begin(), consistent.end(), [](bool c) { return c; })) { continue; }
    if (std::all_of(consistent.begin(), consistent.end(), [](bool c) { return c; })) { continue; }
    for (Size i = 0; i < ids.size(); ++i)
    {
      if (!consistent[i]) { moves.push_back({{ids[i].run->getUuid(), ids[i].query->getId()}, std::nullopt, std::nullopt}); }
    }
  }
  // as unassigned identifications, after the others
  moveIdentifications(multiplets, moves);
  return moves.size();
}
} // namespace OpenMS
