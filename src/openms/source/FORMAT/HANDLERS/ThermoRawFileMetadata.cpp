// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ThermoRawFileMetadata.h>

#include <cmath>
#include <locale>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace OpenMS::Internal
{
std::string ThermoRawFileMetadata::trailer(const ThermoScan& scan, const std::string& label)
{
  for (const auto& entry : scan.trailer)
  {
    if (entry.label == label) { return entry.value; }
  }
  return "";
}

std::optional<double> ThermoRawFileMetadata::number(const std::string& value)
{
  std::istringstream stream(value);
  stream.imbue(std::locale::classic());
  double result;
  if (! (stream >> result) || ! std::isfinite(result)) { return std::nullopt; }
  stream >> std::ws;
  if (! stream.eof()) { return std::nullopt; }
  return result;
}

double ThermoRawFileMetadata::selectedIon(double target, std::optional<double> mono, double width)
{
  if (! mono || *mono <= 1e-8) { return target; }
  const double mz = *mono;
  if (std::abs(mz - target) <= 0.0001) { return target; }
  const double half = width / 2;
  if (half <= 2.0) { return mz >= target - 3.0 && mz <= target + 2.5 ? mz : target; }
  return mz >= target - half && mz <= target + half ? mz : target;
}

std::vector<ThermoPrecursor> ThermoRawFileMetadata::precursors(const ThermoScan& scan)
{
  const int scan_number = scan.scan_number;
  const int level = scan.ms_level;
  const auto& reactions = scan.reactions;
  State state;
  state.level = level;
  state.native_id = scan.native_id;
  if (level == 1)
  {
    filters_[""] = scan_number;
    scans_[scan_number] = state;
    return {};
  }
  std::smatch match;
  std::string key;
  static const std::regex filter_pattern(R"(ms\d+ (.+?) \[)");
  if (std::regex_search(scan.filter, match, filter_pattern)) { key = match[1]; }
  const int fallback = parentFromFilter_(key);
  const auto master = number(trailer(scan, "Master Scan Number:"));
  int parent = master && *master > 0 && *master < scan_number ? static_cast<int>(*master) : fallback;
  auto parent_it = scans_.find(parent);
  // A master scan of the same MS order is a sibling (Tribrid decision trees, or several
  // activations of one isolation such as HCD / ETD / EThcD scans of the same precursor).
  // Its reactions do not apply here; descend from the scan the sibling was triggered by.
  if (parent_it != scans_.end() && parent_it->second.level == level)
  {
    parent = parent_it->second.parent;
    parent_it = scans_.find(parent);
  }
  std::size_t index = parent_it == scans_.end() ? lastReaction_(reactions) : parent_it->second.reactions;
  if (index >= reactions.size() && parent_it != scans_.end())
  {
    parent = fallback;
    parent_it = scans_.find(parent);
    index = parent_it == scans_.end() ? lastReaction_(reactions) : parent_it->second.reactions;
  }
  state.parent = parent_it == scans_.end() ? 0 : parent;
  std::vector<ThermoPrecursor> result;
  if (index < reactions.size())
  {
    const auto& reaction = reactions[index];
    if (! reaction.precursor_mass)
    {
      throw std::runtime_error("Thermo scan " + std::to_string(scan_number) + ": reaction without precursor mass");
    }
    auto width = number(trailer(scan, "MS" + std::to_string(level) + " Isolation Width:"));
    if (! width) { width = reaction.isolation_width; }
    if (width && *width < 0) { width.reset(); }
    const double target = *reaction.precursor_mass;
    const double w = width.value_or(0.0);
    const double offset = reaction.isolation_offset.value_or(0.0);
    const auto charge = number(trailer(scan, "Charge State:"));
    ThermoPrecursor precursor;
    precursor.target_mz = target;
    precursor.selected_mz = selectedIon(target, number(trailer(scan, "Monoisotopic M/Z:")), w);
    if (charge && *charge > 0) { precursor.charge = static_cast<int>(*charge); }
    precursor.width = width;
    precursor.lower_offset = w / 2 - offset;
    precursor.upper_offset = w / 2 + offset;
    precursor.parent_scan = state.parent;
    precursor.estimate_intensity = true;
    precursor.spectrum_ref = parent_it == scans_.end() ? "" : parent_it->second.native_id;
    precursor.activation = reaction;
    ++index;
    // The parent's consumed reactions already account for earlier MS levels.
    // Like TRFP, retain the remaining reaction as supplemental activation even
    // when the vendor's activation flag or type is unexpected.
    if (index < reactions.size()) { precursor.supplemental = reactions[index++]; }
    state.own.push_back(precursor);
    const auto masses = spsMasses_(scan);
    // First SPS selection is represented by the primary reaction.
    for (std::size_t i = 1; i < masses.size(); ++i)
    {
      auto sps = precursor;
      sps.target_mz = masses[i];
      sps.selected_mz = masses[i];
      sps.charge.reset();
      sps.estimate_intensity = false;
      state.own.push_back(std::move(sps));
    }
    state.reactions = index;
    result = state.own;
    // Ancestors are resolved through the parent chain (MS3 -> MS2 -> MS1) instead of being
    // copied into every scan, so retained metadata stays one descriptor per MSn scan.
    for (int ancestor = state.parent; ancestor != 0;)
    {
      const auto it = scans_.find(ancestor);
      if (it == scans_.end()) { break; }
      result.insert(result.end(), it->second.own.begin(), it->second.own.end());
      ancestor = it->second.parent;
    }
  }
  if (! key.empty()) { filters_[key] = scan_number; }
  scans_[scan_number] = std::move(state);
  return result;
}

bool ThermoRawFileMetadata::supplemental_(const ThermoReaction& first, const ThermoReaction& second)
{
  const std::string& a = first.activation;
  const std::string& b = second.activation;
  return first.precursor_mass && second.precursor_mass && std::abs(*first.precursor_mass - *second.precursor_mass) < 0.0001
         && (a == "ElectronTransferDissociation" || a == "ElectronCaptureDissociation")
         && (b == "HigherEnergyCollisionalDissociation" || b == "CollisionInducedDissociation");
}

std::size_t ThermoRawFileMetadata::lastReaction_(const std::vector<ThermoReaction>& reactions)
{
  if (reactions.empty()) { return 0; }
  std::size_t index = reactions.size() - 1;
  if (index > 0 && supplemental_(reactions[index - 1], reactions[index])) { --index; }
  return index;
}

int ThermoRawFileMetadata::parentFromFilter_(const std::string& key) const
{
  std::istringstream stream(key);
  std::vector<std::string> parts;
  for (std::string part; stream >> part;)
  {
    parts.push_back(part);
  }
  if (! parts.empty())
  {
    const auto mass = parts.back().substr(0, parts.back().find('@'));
    while (! parts.empty() && parts.back().substr(0, parts.back().find('@')) == mass)
    {
      parts.pop_back();
    }
  }
  std::string parent_key;
  for (const auto& part : parts)
  {
    if (! parent_key.empty()) { parent_key += ' '; }
    parent_key += part;
  }
  const auto it = filters_.find(parent_key);
  return it == filters_.end() ? 0 : it->second;
}

std::vector<double> ThermoRawFileMetadata::spsMasses_(const ThermoScan& scan)
{
  std::vector<double> result;
  const bool modern = ! trailer(scan, "SPS Masses:").empty();
  static const std::regex legacy_pattern(R"(SPS Mass\s+\d+:)");
  static const std::regex modern_pattern(R"(SPS Masses(?:\s+Continued)?:)");
  for (const auto& entry : scan.trailer)
  {
    if ((! modern && std::regex_match(entry.label, legacy_pattern)) || (modern && std::regex_match(entry.label, modern_pattern)))
    {
      std::istringstream values(entry.value);
      for (std::string value; std::getline(values, value, ',');)
      {
        const auto mass = number(value);
        if (mass && *mass > 0) { result.push_back(*mass); }
      }
    }
  }
  return result;
}
} // namespace OpenMS::Internal
