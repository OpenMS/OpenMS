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
    for (const auto& id : map.getUnassignedPeptideIdentifications())
    {
      cvs.insert(getCV_(id));
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
      part.getUnassignedPeptideIdentifications().clear();
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
    for (const auto& id : maps[i].getUnassignedPeptideIdentifications())
    {
      result.at(getCV_(id))[i].getUnassignedPeptideIdentifications().push_back(id);
    }
  }
  return result;
}

void MS1LabeledFAIMS::append_(ConsensusMap& result, ConsensusMap&& part)
{
  if (result.getColumnHeaders().empty())
  {
    result = std::move(part);
    return;
  }
  // appendRows rewrites filenames and duplicates protein runs. These partitions share columns
  // and protein runs: only feature rows, column sizes and unassigned IDs need concatenating.
  for (const auto& [index, header] : part.getColumnHeaders())
  {
    result.getColumnHeaders().at(index).size += header.size;
  }
  for (auto& feature : part)
  {
    result.push_back(std::move(feature));
  }
  for (auto& id : part.getUnassignedPeptideIdentifications())
  {
    result.getUnassignedPeptideIdentifications().push_back(std::move(id));
  }
}

void MS1LabeledFAIMS::annotate(IDMapper& mapper,
                               ConsensusMap& map,
                               const PeptideIdentificationList& ids,
                               const std::vector<ProteinIdentification>& proteins)
{
  map.setUnassignedPeptideIdentifications(ids);
  auto partitions = split_({map});
  map.clear(true);
  for (auto& [cv, parts] : partitions)
  {
    auto& part = parts.front();
    auto part_ids = std::move(part.getUnassignedPeptideIdentifications());
    part.getUnassignedPeptideIdentifications().clear();
    mapper.annotate(part, part_ids, proteins, true, true);
    append_(map, std::move(part));
  }
  map.sortByPosition();
  map.updateRanges();
}

void MS1LabeledFAIMS::group(FeatureGroupingAlgorithmQT& linker, const std::vector<ConsensusMap>& maps, ConsensusMap& result)
{
  auto partitions = split_(maps);
  result.clear(true);
  for (auto& [cv, parts] : partitions)
  {
    for (auto& part : parts)
    {
      part.updateRanges();
    }
    ConsensusMap linked;
    linker.group(parts, linked);
    linker.transferSubelements(parts, linked);
    if (cv)
    {
      for (auto& feature : linked)
      {
        feature.setMetaValue(Constants::UserParam::FAIMS_CV, *cv);
      }
    }
    append_(result, std::move(linked));
  }
  result.updateRanges();
}
} // namespace OpenMS
