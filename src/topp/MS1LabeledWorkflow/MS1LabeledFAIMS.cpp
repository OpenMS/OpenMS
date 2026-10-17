// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledFAIMS.h"

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmQT.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/METADATA/SpectrumLookup.h>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <set>
#include <unordered_map>
#include <utility>

namespace OpenMS
{
MS1LabeledFAIMS::CV MS1LabeledFAIMS::getCV_(const MetaInfoInterface& value)
{
  if (! value.metaValueExists(Constants::UserParam::FAIMS_CV)) { return std::nullopt; }
  const double cv = value.getMetaValue(Constants::UserParam::FAIMS_CV);
  if (! std::isfinite(cv))
  {
    throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "FAIMS compensation voltages must be finite.");
  }
  return cv;
}

void MS1LabeledFAIMS::annotateCompensationVoltages(const MSExperiment& spectra, PeptideIdentificationList& ids)
{
  std::set<double> ms1_cvs;
  std::vector<CV> spectrum_cvs(spectra.size());
  CV previous_cv;
  for (Size i = 0; i < spectra.size(); ++i)
  {
    const auto& spectrum = spectra[i];
    if (spectrum.getDriftTimeUnit() == DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
    {
      const double cv = spectrum.getDriftTime();
      if (! std::isfinite(cv) || cv == IMTypes::DRIFTTIME_NOT_SET)
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "A FAIMS spectrum has no valid compensation voltage.");
      }
      previous_cv = cv;
      spectrum_cvs[i] = cv;
      if (spectrum.getMSLevel() == 1) { ms1_cvs.insert(cv); }
    }
    else if (spectrum.getMSLevel() > 1)
    {
      // Same acquisition-order fallback as IMDataConverter::splitByFAIMSCV.
      spectrum_cvs[i] = previous_cv;
    }
    else
    {
      previous_cv.reset();
    }
  }
  if (ms1_cvs.empty()) { return; }

  // Missing references may be repaired by RT only after the CV is known. A global RT lookup
  // can pick another CV and collapse distinct spectra during subsequent PSM deduplication.
  std::map<double, std::map<double, std::string>> ms2_references;
  if (std::any_of(ids.begin(), ids.end(), [](const auto& id) { return id.getSpectrumReference().empty(); }))
  {
    for (Size i = 0; i < spectra.size(); ++i)
    {
      if (spectra[i].getMSLevel() != 2 || ! spectrum_cvs[i]) { continue; }
      auto& references = ms2_references[*spectrum_cvs[i]];
      const auto [it, inserted] = references.emplace(spectra[i].getRT(), spectra[i].getNativeID());
      if (! inserted) { it->second.clear(); } // simultaneous spectra at the same CV are ambiguous
    }
  }

  SpectrumLookup lookup;
  lookup.readSpectra(spectra.getSpectra());
  lookup.addReferenceFormat("^scan=(?<SCAN>\\d+)$");
  lookup.addReferenceFormat("^index=(?<INDEX0>\\d+)$");
  for (auto& id : ids)
  {
    CV cv = getCV_(id);
    CV spectrum_cv;
    const auto& reference = id.getSpectrumReference();
    if (! reference.empty())
    {
      std::optional<Size> index;
      try
      {
        index = lookup.findByNativeID(reference);
      }
      catch (const Exception::ElementNotFound&)
      {
        try
        {
          index = lookup.findByReference(reference);
        }
        catch (const Exception::ElementNotFound&)
        {
        }
        catch (const Exception::ParseError&)
        {
        }
      }
      if (index) { spectrum_cv = spectrum_cvs[*index]; }
    }
    if (spectrum_cv)
    {
      if (cv && std::abs(*cv - *spectrum_cv) > 0.01)
      {
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Identification and spectrum disagree on the FAIMS compensation voltage: '" + reference + "'.");
      }
      cv = spectrum_cv;
    }
    // A single-CV acquisition is unambiguous even when references were not exported by the search.
    if (! cv && ms1_cvs.size() == 1) { cv = *ms1_cvs.begin(); }
    if (! cv)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Cannot determine the FAIMS compensation voltage for identification '" + reference
                                            + "'. Multi-CV data require valid spectrum references or FAIMS_CV annotations; RT alone is ambiguous.");
    }
    // Canonicalize independently rounded ID metadata to the actual acquisition CV.
    auto matching_cv = ms1_cvs.end();
    for (auto it = ms1_cvs.begin(); it != ms1_cvs.end(); ++it)
    {
      if (std::abs(*it - *cv) <= 0.01)
      {
        if (matching_cv != ms1_cvs.end())
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Ambiguous FAIMS compensation voltage annotation.");
        }
        matching_cv = it;
      }
    }
    if (matching_cv == ms1_cvs.end())
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Identification '" + reference + "' has a FAIMS compensation voltage absent from the MS1 spectra.");
    }
    id.setMetaValue(Constants::UserParam::FAIMS_CV, *matching_cv);
    if (reference.empty())
    {
      const auto& references = ms2_references[*matching_cv];
      auto match = references.lower_bound(id.getRT());
      if (match != references.begin() && (match == references.end() || id.getRT() - std::prev(match)->first < match->first - id.getRT())) { --match; }
      if (match == references.end() || ! std::isfinite(id.getRT()) || std::abs(match->first - id.getRT()) > lookup.rt_tolerance
          || match->second.empty())
      {
        throw Exception::MissingInformation(
          __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Cannot repair a missing spectrum reference at the identification's FAIMS CV and RT. Provide valid MS2 spectrum references.");
      }
      id.setSpectrumReference(match->second);
    }
  }
}

std::map<MS1LabeledFAIMS::CV, std::vector<ConsensusMap>> MS1LabeledFAIMS::split_(const std::vector<ConsensusMap>& maps)
{
  std::set<CV> cvs;
  for (const auto& map : maps)
  {
    for (const auto& feature : map)
    {
      cvs.insert(getCV_(feature));
    }
  }
  if (cvs.empty()) { cvs.insert(std::nullopt); }

  std::map<CV, std::vector<ConsensusMap>> result;
  for (const auto cv : cvs)
  {
    auto& parts = result[cv];
    parts = maps;
    for (auto& part : parts)
    {
      part.resize(0);
      for (auto& [index, header] : part.getColumnHeaders())
      {
        header.size = 0;
      }
    }
  }
  for (Size i = 0; i < maps.size(); ++i)
  {
    for (const auto& feature : maps[i])
    {
      auto& part = result.at(getCV_(feature))[i];
      part.push_back(feature);
      for (const auto& handle : feature.getFeatures())
      {
        ++part.getColumnHeaders().at(handle.getMapIndex()).size;
      }
    }
  }
  return result;
}

void MS1LabeledFAIMS::annotate(const IDMapper& mapper, ConsensusMap& map, const FeatureMap& channels, const IdentificationData& ids)
{
  using ID = IdentificationData;

  // channel feature (= handle unique id) -> (multiplet, channel)
  std::unordered_map<UInt64, std::pair<Size, Size>> channel_of;
  for (Size m = 0; m < map.size(); ++m)
  {
    for (const FeatureHandle& handle : map[m].getFeatures())
    {
      channel_of.emplace(handle.getUniqueId(), std::make_pair(m, static_cast<Size>(handle.getMapIndex())));
    }
  }

  // The channels and IDs of one CV are mapped together; the channels get the identification data of the map, which
  // the mapping adds the runs of the IDs to.
  IDMapper constrained = mapper;
  Param param = constrained.getParameters();
  param.setValue("match_meta_value", Constants::UserParam::FAIMS_CV);
  constrained.setParameters(param);
  FeatureMap part = channels;
  part.getIdentificationData() = std::move(map.getIdentificationData());
  part.updateRanges();
  const Size recorded = part.getDataProcessing().size();
  // RT and m/z from the convex hulls: the elution region of the channel, plus the tolerances
  constrained.annotate(part, ids, false, false);
  map.getDataProcessing().insert(map.getDataProcessing().end(), part.getDataProcessing().begin() + recorded, part.getDataProcessing().end());
  ID& data = part.getIdentificationData();

  // the source of every identification, for its copies
  std::map<ID::QueryReference, ID::SourceId> source_of;
  for (const auto& run : data.getRuns())
  {
    for (const auto& source : run.getSources())
    {
      for (const auto& query : source.identifications) { source_of.emplace(ID::QueryReference {run.getUuid(), query.getId()}, source.id); }
    }
  }

  // Every assignment of an ID to a multiplet is a copy of it, with the matched channel as map index. A channel links
  // the IDs it matched with all their matches.
  struct Copy
  {
    std::string run;
    ID::SourceId source;
    ID::Observation observation;
    std::vector<std::pair<ID::MatchData, std::vector<std::optional<double>>>> matches;
    Size multiplet;
  };
  std::vector<Copy> copies;
  std::set<ID::QueryReference> copied;
  std::vector<std::set<ID::QueryReference>> attached(map.size()); // per multiplet: the IDs attached to it
  for (const Feature& channel : part)
  {
    const auto linked = channel.getLinkedIdentifications(data);
    if (linked.empty()) { continue; }
    const auto it = channel_of.find(channel.getUniqueId());
    if (it == channel_of.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "A channel feature does not belong to any multiplet.",
                                    StringUtils::toStr(channel.getUniqueId()));
    }
    const auto [multiplet, map_index] = it->second;
    for (const auto& entry : linked)
    {
      const ID::QueryReference reference {entry.run->getUuid(), entry.query->getId()};
      if (!attached[multiplet].insert(reference).second) { continue; } // already attached through another channel
      Copy copy {entry.run->getIdentifier(), source_of.at(reference), entry.query->getObservation(), {}, multiplet};
      copy.observation.setMetaValue("map_index", map_index);
      for (const ID::Match* match : entry.matches) { copy.matches.emplace_back(match->getData(), entry.run->getScores(*match)); }
      copies.push_back(std::move(copy));
      copied.insert(reference);
    }
  }

  // the copies are listed after all other identifications, in the order of assignment
  UInt64 next = 1;
  for (const auto& run : data.getRuns()) { next = std::max(next, run.getNextQueryId()); }
  for (auto& copy : copies)
  {
    auto& run = data.getRun(copy.run);
    const ID::QueryId query = run.importIdentification(copy.source, ID::QueryId {next++}, std::move(copy.observation));
    ConsensusFeature& feature = map[copy.multiplet];
    feature.addIDQuery({run.getUuid(), query});
    for (const auto& [match, scores] : copy.matches) { feature.addIDMatch({run.getUuid(), run.addMatch(query, match, scores)}); }
  }
  // the attached IDs are replaced by their copies; the others stay unassigned
  data.eraseIdentifications([&](const ID::Run& run, const ID::Identification& query) { return copied.contains({run.getUuid(), query.getId()}); });
  map.getIdentificationData() = std::move(data);
  map.updateRanges();
}

void MS1LabeledFAIMS::group(FeatureGroupingAlgorithmQT& linker, const std::vector<ConsensusMap>& maps, ConsensusMap& result)
{
  auto partitions = split_(maps);
  result.clear(true);
  bool first = true;
  for (auto& [cv, parts] : partitions)
  {
    for (auto& part : parts)
    {
      part.updateRanges();
    }
    ConsensusMap linked;
    linker.group(parts, linked);
    if (cv)
    {
      for (auto& feature : linked)
      {
        feature.setMetaValue(Constants::UserParam::FAIMS_CV, *cv);
      }
    }
    // The partitions share protein runs; the columns follow from the input maps below.
    if (first)
    {
      result = std::move(linked);
      first = false;
      continue;
    }
    for (auto& feature : linked)
    {
      result.push_back(std::move(feature));
    }
  }
  // Every partition has all identifications of its maps (only the features are split), and linking a partition drops
  // those that only its features left out link. The identifications come from the input maps instead, for the features
  // of all partitions.
  if (std::any_of(maps.begin(), maps.end(), [](const ConsensusMap& map) { return ! map.getIdentificationData().empty(); }))
  {
    FeatureGroupingAlgorithm::groupIdentifications(maps, result);
  }
  // the partitions have the features and columns of the input maps: their channels are the subelements of those
  linker.transferSubelements(maps, result);
  // as in the partitions, a column has the size of its channel values
  Size column = 0;
  for (const auto& map : maps)
  {
    std::map<UInt64, Size> sizes;
    for (const auto& feature : map)
    {
      for (const auto& handle : feature.getFeatures())
      {
        ++sizes[handle.getMapIndex()];
      }
    }
    for (const auto& [index, header] : map.getColumnHeaders())
    {
      result.getColumnHeaders().at(column++).size = sizes[index];
    }
  }
  result.updateRanges();
}
} // namespace OpenMS
