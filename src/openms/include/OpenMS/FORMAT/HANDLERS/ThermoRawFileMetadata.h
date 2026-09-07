// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------
#pragma once

#include <cmath>
#include <locale>
#include <map>
#include <nlohmann/json.hpp>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace OpenMS::Internal
{
/** @brief Resolve the precursor hierarchy from bridge metadata without fetching
   peaks. Keeps isolation targets distinct from selected ions, and retains SPS
   selections and supplemental reactions. The filter fallback follows
   ThermoRawFileParser. */
class ThermoRawFileMetadata
{
public:
  static std::string trailer(const nlohmann::json& scan, const std::string& label)
  {
    for (const auto& entry : scan.at("trailer"))
    {
      if (entry.at("label") == label && entry.at("value").is_string()) { return entry.at("value").get<std::string>(); }
    }
    return "";
  }

  static nlohmann::json number(const std::string& value)
  {
    std::istringstream stream(value);
    stream.imbue(std::locale::classic());
    double result;
    if (! (stream >> result) || ! std::isfinite(result)) { return nullptr; }
    stream >> std::ws;
    return stream.eof() ? nlohmann::json(result) : nlohmann::json(nullptr);
  }

  static double selectedIon(double target, const nlohmann::json& mono, double width)
  {
    if (! mono.is_number() || mono.get<double>() <= 1e-8) { return target; }
    double mz = mono.get<double>();
    if (std::abs(mz - target) <= 0.0001) { return target; }
    double half = width / 2;
    if (half <= 2.0) { return mz >= target - 3.0 && mz <= target + 2.5 ? mz : target; }
    return mz >= target - half && mz <= target + half ? mz : target;
  }

  nlohmann::json precursors(const nlohmann::json& scan)
  {
    const int scan_number = scan.at("scan_number").get<int>();
    const int level = scan.at("ms_level").get<int>();
    const auto& reactions = scan.at("reactions");
    State state;
    state.level = level;
    state.native_id = scan.at("native_id").get<std::string>();
    if (level == 1)
    {
      filters_[""] = scan_number;
      scans_[scan_number] = state;
      return state.precursors;
    }
    std::smatch match;
    const std::string filter = scan.at("filter").get<std::string>();
    std::string key;
    static const std::regex filter_pattern(R"(ms\d+ (.+?) \[)");
    if (std::regex_search(filter, match, filter_pattern)) { key = match[1]; }
    int fallback = parentFromFilter_(key);
    auto master = number(trailer(scan, "Master Scan Number:"));
    int parent
      = master.is_number() && master.get<double>() > 0 && master.get<double>() < scan_number ? static_cast<int>(master.get<double>()) : fallback;
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
    if (index < reactions.size())
    {
      const auto& reaction = reactions[index];
      auto width = number(trailer(scan, "MS" + std::to_string(level) + " Isolation Width:"));
      if (width.is_null()) { width = reaction.at("isolation_width"); }
      if (width.is_number() && width.get<double>() < 0) { width = nullptr; }
      const double target = reaction.at("precursor_mass").get<double>();
      const double w = width.is_number() ? width.get<double>() : 0.0;
      const double offset = reaction.at("isolation_offset").is_number() ? reaction.at("isolation_offset").get<double>() : 0.0;
      auto charge = number(trailer(scan, "Charge State:"));
      if (charge.is_number() && charge.get<double>() <= 0) { charge = nullptr; }
      nlohmann::json precursor = {{"target_mz", target},
                                  {"selected_mz", selectedIon(target, number(trailer(scan, "Monoisotopic M/Z:")), w)},
                                  {"charge", charge},
                                  {"width", width},
                                  {"lower_offset", w / 2 - offset},
                                  {"upper_offset", w / 2 + offset},
                                  {"parent_scan", parent_it == scans_.end() ? 0 : parent},
                                  {"estimate_intensity", true},
                                  {"spectrum_ref", parent_it == scans_.end() ? "" : parent_it->second.native_id},
                                  {"activation", reaction},
                                  {"supplemental", nullptr}};
      ++index;
      // The parent's consumed reactions already account for earlier MS levels.
      // Like TRFP, retain the remaining reaction as supplemental activation even
      // when the vendor's activation flag or type is unexpected.
      if (index < reactions.size()) { precursor["supplemental"] = reactions[index++]; }
      state.precursors.push_back(precursor);
      auto masses = spsMasses_(scan);
      // First SPS selection is represented by the primary reaction.
      for (std::size_t i = 1; i < masses.size(); ++i)
      {
        auto sps = precursor;
        sps["target_mz"] = masses[i];
        sps["selected_mz"] = masses[i];
        sps["charge"] = nullptr;
        sps["estimate_intensity"] = false;
        state.precursors.push_back(std::move(sps));
      }
      if (parent_it != scans_.end())
      {
        for (const auto& ancestor : parent_it->second.precursors)
        {
          state.precursors.push_back(ancestor);
        }
      }
      state.reactions = index;
    }
    if (! key.empty()) { filters_[key] = scan_number; }
    scans_[scan_number] = state;
    return state.precursors;
  }

private:
  struct State
  {
    int level = 1;
    int parent = 0;              ///< scan this one descends from (0 if unknown)
    std::size_t reactions = 0;   ///< reactions consumed up to and including this scan
    std::string native_id;
    nlohmann::json precursors = nlohmann::json::array();
  };
  std::map<int, State> scans_;
  std::map<std::string, int> filters_;

  static bool supplemental_(const nlohmann::json& first, const nlohmann::json& second)
  {
    const std::string a = first.at("activation").get<std::string>(), b = second.at("activation").get<std::string>();
    return first.at("precursor_mass").is_number() && second.at("precursor_mass").is_number()
           && std::abs(first.at("precursor_mass").get<double>() - second.at("precursor_mass").get<double>()) < 0.0001
           && (a == "ElectronTransferDissociation" || a == "ElectronCaptureDissociation")
           && (b == "HigherEnergyCollisionalDissociation" || b == "CollisionInducedDissociation");
  }
  static std::size_t lastReaction_(const nlohmann::json& reactions)
  {
    if (reactions.empty()) { return 0; }
    std::size_t index = reactions.size() - 1;
    if (index > 0 && supplemental_(reactions[index - 1], reactions[index])) { --index; }
    return index;
  }
  int parentFromFilter_(const std::string& key) const
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
    auto it = filters_.find(parent_key);
    return it == filters_.end() ? 0 : it->second;
  }
  static std::vector<double> spsMasses_(const nlohmann::json& scan)
  {
    std::vector<double> result;
    const bool modern = ! trailer(scan, "SPS Masses:").empty();
    static const std::regex legacy_pattern(R"(SPS Mass\s+\d+:)");
    static const std::regex modern_pattern(R"(SPS Masses(?:\s+Continued)?:)");
    for (const auto& entry : scan.at("trailer"))
    {
      const std::string label = entry.at("label").get<std::string>();
      if ((! modern && std::regex_match(label, legacy_pattern)) || (modern && std::regex_match(label, modern_pattern)))
      {
        std::istringstream values(entry.at("value").get<std::string>());
        for (std::string value; std::getline(values, value, ',');)
        {
          auto mass = number(value);
          if (mass.is_number() && mass.get<double>() > 0) { result.push_back(mass.get<double>()); }
        }
      }
    }
    return result;
  }
};
} // namespace OpenMS::Internal
