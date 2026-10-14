// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/ThermoRawFileMetadata.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>

#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>

namespace OpenMS::Internal
{
namespace
{
  /// Whitespace as the classic locale defines it (what istream >> skips).
  bool isSpace(char c)
  {
    return c == ' ' || c == '\t' || c == '\n' || c == '\v' || c == '\f' || c == '\r';
  }

  bool isDigit(char c) { return c >= '0' && c <= '9'; }
} // namespace

ThermoTrailerIndex::Key ThermoTrailerIndex::add(const std::string& label)
{
  labels_.push_back(label);
  positions_.push_back(npos_);
  return labels_.size() - 1;
}

const std::string* ThermoTrailerIndex::find(const ThermoScan& scan, Key key)
{
  const std::string& label = labels_[key];
  std::size_t& position = positions_[key];
  if (position < scan.trailer.size() && scan.trailer[position].label == label) { return &scan.trailer[position].value; }
  for (std::size_t i = 0; i < scan.trailer.size(); ++i)
  {
    if (scan.trailer[i].label == label)
    {
      position = i;
      return &scan.trailer[i].value;
    }
  }
  return nullptr;
}

const std::string& ThermoTrailerIndex::value(const ThermoScan& scan, Key key)
{
  static const std::string empty;
  const std::string* result = find(scan, key);
  return result ? *result : empty;
}

ThermoRawFileMetadata::ThermoRawFileMetadata() :
  master_key_(trailer_.add("Master Scan Number:")),
  charge_key_(trailer_.add("Charge State:")),
  mono_key_(trailer_.add("Monoisotopic M/Z:")),
  sps_masses_key_(trailer_.add("SPS Masses:"))
{
}

std::string ThermoRawFileMetadata::trailer(const ThermoScan& scan, const std::string& label)
{
  for (const auto& entry : scan.trailer)
  {
    if (entry.label == label) { return entry.value; }
  }
  return "";
}

std::optional<double> ThermoRawFileMetadata::number(std::string_view value)
{
  const char* begin = value.data();
  const char* end = begin + value.size();
  while (begin != end && isSpace(*begin)) { ++begin; }
  while (end != begin && isSpace(*(end - 1))) { --end; }
  // A single leading '+' is allowed (as by istream), but not a second sign.
  if (begin != end && *begin == '+' && (begin + 1 == end || begin[1] == '+' || begin[1] == '-')) { return std::nullopt; }
  double result;
  if (begin == end || ! StringUtilsHelper::extractDouble(begin, end, result) || begin != end || ! std::isfinite(result))
  {
    return std::nullopt;
  }
  return result;
}

std::string ThermoRawFileMetadata::filterKey(const std::string& filter)
{
  // Equivalent of std::regex_search with R"(ms\d+ (.+?) \[)" ('.' excludes line terminators).
  for (std::size_t p = filter.find("ms"); p != std::string::npos; p = filter.find("ms", p + 1))
  {
    std::size_t q = p + 2;
    if (q >= filter.size() || ! isDigit(filter[q])) { continue; }
    while (q < filter.size() && isDigit(filter[q])) { ++q; }
    if (q >= filter.size() || filter[q] != ' ') { continue; }
    const std::size_t start = q + 1;
    if (start >= filter.size()) { continue; }
    const std::size_t stop = filter.find(" [", start + 1);
    if (stop == std::string::npos) { return ""; } // no later start can match either
    if (filter.find_first_of("\r\n", start) < stop) { continue; }
    return filter.substr(start, stop - start);
  }
  return "";
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
  if (first_seen_ == 0) { first_seen_ = scan_number; }
  prune_(scan_number);
  State state;
  state.level = level;
  state.native_id = scan.native_id;
  if (level == 1)
  {
    filters_[""] = scan_number;
    scans_[scan_number] = state;
    return {};
  }
  const std::string key = filterKey(scan.filter);
  const int fallback = parentFromFilter_(key);
  const auto master = number(trailer_.value(scan, master_key_));
  int parent = master && *master > 0 && *master < scan_number ? static_cast<int>(*master) : fallback;
  if (parent == static_cast<int>(master.value_or(0)) && parent >= first_seen_ && parent < forgotten_below_ && ! warned_forgotten_
      && scans_.find(parent) == scans_.end())
  {
    warned_forgotten_ = true;
    OPENMS_LOG_WARN << "Thermo RAW: scan " << scan_number << " references master scan " << parent << ", more than "
                    << retentionWindow() << " scans before it; its precursor is resolved without the master scan.\n";
  }
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
    if (width_keys_.size() <= static_cast<std::size_t>(level))
    {
      for (std::size_t l = width_keys_.size(); l <= static_cast<std::size_t>(level); ++l)
      {
        width_keys_.push_back(trailer_.add("MS" + std::to_string(l) + " Isolation Width:"));
      }
    }
    auto width = number(trailer_.value(scan, width_keys_[level]));
    if (! width) { width = reaction.isolation_width; }
    if (width && *width < 0) { width.reset(); }
    const double target = *reaction.precursor_mass;
    const double w = width.value_or(0.0);
    const double offset = reaction.isolation_offset.value_or(0.0);
    const auto charge = number(trailer_.value(scan, charge_key_));
    ThermoPrecursor precursor;
    precursor.target_mz = target;
    precursor.selected_mz = selectedIon(target, number(trailer_.value(scan, mono_key_)), w);
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

void ThermoRawFileMetadata::prune_(int current)
{
  // Prune in batches, so the amortized cost per scan stays constant.
  if (scans_.size() < 2 * static_cast<std::size_t>(retentionWindow())) { return; }
  const int cutoff = current - retentionWindow();
  std::set<int> keep;
  auto keepChain = [&](int scan) {
    while (scan != 0 && keep.insert(scan).second)
    {
      const auto it = scans_.find(scan);
      if (it == scans_.end()) { break; }
      scan = it->second.parent;
    }
  };
  for (auto it = filters_.begin(); it != filters_.end();)
  {
    // The most recent MS1 scan stays known, however long ago it was acquired.
    if (it->second < cutoff && ! it->first.empty()) { it = filters_.erase(it); }
    else
    {
      keepChain(it->second);
      ++it;
    }
  }
  for (auto it = scans_.lower_bound(cutoff); it != scans_.end(); ++it)
  {
    keepChain(it->first);
  }
  for (auto it = scans_.begin(); it != scans_.end() && it->first < cutoff;)
  {
    if (keep.count(it->first) == 0) { it = scans_.erase(it); }
    else { ++it; }
  }
  forgotten_below_ = cutoff;
}

int ThermoRawFileMetadata::parentFromFilter_(const std::string& key) const
{
  std::vector<std::string_view> parts;
  for (std::size_t i = 0; i < key.size();)
  {
    while (i < key.size() && isSpace(key[i])) { ++i; }
    const std::size_t start = i;
    while (i < key.size() && ! isSpace(key[i])) { ++i; }
    if (i > start) { parts.emplace_back(key.data() + start, i - start); }
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

bool ThermoRawFileMetadata::isLegacySpsLabel_(const std::string& label)
{
  // R"(SPS Mass\s+\d+:)"
  static constexpr std::string_view prefix = "SPS Mass";
  if (label.compare(0, prefix.size(), prefix) != 0) { return false; }
  std::size_t i = prefix.size();
  const std::size_t spaces = i;
  while (i < label.size() && isSpace(label[i])) { ++i; }
  if (i == spaces) { return false; }
  const std::size_t digits = i;
  while (i < label.size() && isDigit(label[i])) { ++i; }
  return i > digits && i + 1 == label.size() && label[i] == ':';
}

bool ThermoRawFileMetadata::isModernSpsLabel_(const std::string& label)
{
  // R"(SPS Masses(?:\s+Continued)?:)"
  static constexpr std::string_view prefix = "SPS Masses";
  if (label.compare(0, prefix.size(), prefix) != 0) { return false; }
  std::size_t i = prefix.size();
  if (i + 1 == label.size() && label[i] == ':') { return true; }
  const std::size_t spaces = i;
  while (i < label.size() && isSpace(label[i])) { ++i; }
  return i > spaces && label.compare(i, std::string::npos, "Continued:") == 0;
}

std::vector<double> ThermoRawFileMetadata::spsMasses_(const ThermoScan& scan)
{
  std::vector<double> result;
  const bool modern = ! trailer_.value(scan, sps_masses_key_).empty();
  for (const auto& entry : scan.trailer)
  {
    if (entry.label.size() < 9 || entry.label[0] != 'S') { continue; } // fast reject before the full match
    if ((! modern && isLegacySpsLabel_(entry.label)) || (modern && isModernSpsLabel_(entry.label)))
    {
      const std::string_view values(entry.value);
      for (std::size_t start = 0; start <= values.size();)
      {
        std::size_t stop = values.find(',', start);
        if (stop == std::string_view::npos) { stop = values.size(); }
        const auto mass = number(values.substr(start, stop - start));
        if (mass && *mass > 0) { result.push_back(*mass); }
        start = stop + 1;
      }
    }
  }
  return result;
}
} // namespace OpenMS::Internal
