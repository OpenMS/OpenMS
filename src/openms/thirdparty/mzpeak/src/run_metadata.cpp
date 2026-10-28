/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/run_metadata.h"

#include <cstdint>
#include <limits>

#include "mzpeak/util/json_writer.h"

/*
 * RDR-24 — parse the run-level metadata blocks from the index `metadata{}`
 * object into the typed model declared in run_metadata.h.  Boost.JSON is used
 * throughout to match the `files[]` parsing in src/index.cpp and src/schema.
 *
 * The parsing is intentionally tolerant: every field is optional, missing or
 * mistyped members are skipped rather than throwing, since this metadata is
 * descriptive (not load-bearing for reading the binary arrays) and the spec is
 * a pre-1.0 living standard.
 */
namespace MzPeak {

// Local alias: the public header qualifies Boost.JSON in full; inside this
// translation unit the short form keeps the parsing code readable.
namespace json = boost::json;

/******************************************************************************/
namespace {

/// Read an optional string member; absent or JSON-null yields std::nullopt.
std::optional<std::string> opt_string(const json::object& o, const char* key)
{
  if (const auto it = o.find(key); it != o.end() && it->value().is_string()) {
    return std::string(it->value().as_string().c_str());
  }
  return std::nullopt;
}

/// Read an optional integer member; absent/null/non-integral yields nullopt.
/// Accepts both signed and (in-range) unsigned JSON integers, since Boost.JSON
/// parses a non-negative literal as `uint64` and would otherwise be dropped.
std::optional<std::int64_t> opt_int(const json::object& o, const char* key)
{
  const auto it = o.find(key);
  if (it == o.end()) return std::nullopt;
  const auto& v = it->value();
  if (v.is_int64()) return v.as_int64();
  if (v.is_uint64()) {
    const auto u = v.as_uint64();
    if (u <= static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
      return static_cast<std::int64_t>(u);
    }
  }
  return std::nullopt;
}

/// Stringify a CV-param `value` member, which the index mixes between strings
/// and numbers (e.g. `"SN06061F"` vs `1`).  JSON-null yields std::nullopt.
std::optional<std::string> value_to_string(const json::value& v)
{
  switch (v.kind()) {
  case json::kind::string:
    return std::string(v.as_string().c_str());
  case json::kind::int64:
    return std::to_string(v.as_int64());
  case json::kind::uint64:
    return std::to_string(v.as_uint64());
  case json::kind::double_:
    return std::to_string(v.as_double());
  case json::kind::bool_:
    return v.as_bool() ? std::string("true") : std::string("false");
  default:
    return std::nullopt;
  }
}

/// Read an array member as a JSON array (empty if absent or not an array).
json::array array_member(const json::object& o, const char* key)
{
  if (const auto it = o.find(key); it != o.end() && it->value().is_array()) {
    return it->value().as_array();
  }
  return json::array();
}

/// WRT-2 — set an optional string member, omitting it entirely when absent so
/// `parse(serialize(x)) == x` holds (the parser tolerates missing members).
void put_opt(json::object& o, const char* key, const std::optional<std::string>& v)
{
  if (v.has_value()) o[key] = *v;
}

/// WRT-2 — set an optional integer member, omitting it when absent.
void put_opt(json::object& o, const char* key, const std::optional<std::int64_t>& v)
{
  if (v.has_value()) o[key] = *v;
}

} // namespace

/******************************************************************************/
CvParam CvParam::from_json(const json::object& o)
{
  CvParam p;
  p.accession = opt_string(o, "accession");
  p.name = opt_string(o, "name");
  p.unit = opt_string(o, "unit");
  if (const auto it = o.find("value"); it != o.end()) {
    p.value = value_to_string(it->value());
  }
  return p;
}

/******************************************************************************/
json::object CvParam::to_json() const
{
  json::object o;
  put_opt(o, "accession", accession);
  put_opt(o, "name", name);
  put_opt(o, "unit", unit);
  put_opt(o, "value", value);
  return o;
}

/******************************************************************************/
std::vector<CvParam> cv_params_from_json(const json::array& a)
{
  std::vector<CvParam> params;
  params.reserve(a.size());
  for (const auto& v : a) {
    if (v.is_object()) params.push_back(CvParam::from_json(v.as_object()));
  }
  return params;
}

/******************************************************************************/
json::array cv_params_to_json(const std::vector<CvParam>& params)
{
  json::array a;
  a.reserve(params.size());
  for (const auto& p : params)
    a.push_back(p.to_json());
  return a;
}

/******************************************************************************/
Software Software::from_json(const json::object& o)
{
  Software s;
  if (auto id = opt_string(o, "id")) s.id = std::move(*id);
  s.version = opt_string(o, "version");
  s.parameters = cv_params_from_json(array_member(o, "parameters"));
  return s;
}

/******************************************************************************/
json::object Software::to_json() const
{
  json::object o;
  o["id"] = id;
  put_opt(o, "version", version);
  o["parameters"] = cv_params_to_json(parameters);
  return o;
}

/******************************************************************************/
Component Component::from_json(const json::object& o)
{
  Component c;
  c.component_type = opt_string(o, "component_type");
  c.order = opt_int(o, "order");
  c.parameters = cv_params_from_json(array_member(o, "parameters"));
  return c;
}

/******************************************************************************/
json::object Component::to_json() const
{
  json::object o;
  put_opt(o, "component_type", component_type);
  put_opt(o, "order", order);
  o["parameters"] = cv_params_to_json(parameters);
  return o;
}

/******************************************************************************/
InstrumentConfiguration InstrumentConfiguration::from_json(const json::object& o)
{
  InstrumentConfiguration ic;
  ic.id = opt_int(o, "id");
  ic.software_reference = opt_string(o, "software_reference");
  ic.parameters = cv_params_from_json(array_member(o, "parameters"));
  for (const auto& v : array_member(o, "components")) {
    if (v.is_object()) ic.components.push_back(Component::from_json(v.as_object()));
  }
  return ic;
}

/******************************************************************************/
json::object InstrumentConfiguration::to_json() const
{
  json::object o;
  put_opt(o, "id", id);
  put_opt(o, "software_reference", software_reference);
  o["parameters"] = cv_params_to_json(parameters);
  json::array comps;
  comps.reserve(components.size());
  for (const auto& c : components)
    comps.push_back(c.to_json());
  o["components"] = std::move(comps);
  return o;
}

/******************************************************************************/
ProcessingMethod ProcessingMethod::from_json(const json::object& o)
{
  ProcessingMethod m;
  m.order = opt_int(o, "order");
  m.software_reference = opt_string(o, "software_reference");
  m.parameters = cv_params_from_json(array_member(o, "parameters"));
  return m;
}

/******************************************************************************/
json::object ProcessingMethod::to_json() const
{
  json::object o;
  put_opt(o, "order", order);
  put_opt(o, "software_reference", software_reference);
  o["parameters"] = cv_params_to_json(parameters);
  return o;
}

/******************************************************************************/
DataProcessing DataProcessing::from_json(const json::object& o)
{
  DataProcessing dp;
  if (auto id = opt_string(o, "id")) dp.id = std::move(*id);
  for (const auto& v : array_member(o, "methods")) {
    if (v.is_object())
      dp.methods.push_back(ProcessingMethod::from_json(v.as_object()));
  }
  return dp;
}

/******************************************************************************/
json::object DataProcessing::to_json() const
{
  json::object o;
  o["id"] = id;
  json::array ms;
  ms.reserve(methods.size());
  for (const auto& m : methods)
    ms.push_back(m.to_json());
  o["methods"] = std::move(ms);
  return o;
}

/******************************************************************************/
Sample Sample::from_json(const json::object& o)
{
  Sample s;
  if (auto id = opt_string(o, "id")) s.id = std::move(*id);
  s.name = opt_string(o, "name");
  s.parameters = cv_params_from_json(array_member(o, "parameters"));
  return s;
}

/******************************************************************************/
json::object Sample::to_json() const
{
  json::object o;
  o["id"] = id;
  put_opt(o, "name", name);
  o["parameters"] = cv_params_to_json(parameters);
  return o;
}

/******************************************************************************/
SourceFile SourceFile::from_json(const json::object& o)
{
  SourceFile sf;
  if (auto id = opt_string(o, "id")) sf.id = std::move(*id);
  sf.name = opt_string(o, "name");
  sf.location = opt_string(o, "location");
  sf.parameters = cv_params_from_json(array_member(o, "parameters"));
  return sf;
}

/******************************************************************************/
json::object SourceFile::to_json() const
{
  json::object o;
  o["id"] = id;
  put_opt(o, "name", name);
  put_opt(o, "location", location);
  o["parameters"] = cv_params_to_json(parameters);
  return o;
}

/******************************************************************************/
FileDescription FileDescription::from_json(const json::object& o)
{
  FileDescription fd;
  fd.contents = cv_params_from_json(array_member(o, "contents"));
  for (const auto& v : array_member(o, "source_files")) {
    if (v.is_object())
      fd.source_files.push_back(SourceFile::from_json(v.as_object()));
  }
  return fd;
}

/******************************************************************************/
json::object FileDescription::to_json() const
{
  json::object o;
  o["contents"] = cv_params_to_json(contents);
  json::array sfs;
  sfs.reserve(source_files.size());
  for (const auto& sf : source_files)
    sfs.push_back(sf.to_json());
  o["source_files"] = std::move(sfs);
  return o;
}

/******************************************************************************/
Run Run::from_json(const json::object& o)
{
  Run r;
  r.id = opt_string(o, "id");
  r.start_time = opt_string(o, "start_time");
  r.default_instrument_id = opt_int(o, "default_instrument_id");
  r.default_data_processing_id = opt_string(o, "default_data_processing_id");
  r.default_source_file_id = opt_string(o, "default_source_file_id");
  return r;
}

/******************************************************************************/
json::object Run::to_json() const
{
  json::object o;
  put_opt(o, "id", id);
  put_opt(o, "start_time", start_time);
  put_opt(o, "default_instrument_id", default_instrument_id);
  put_opt(o, "default_data_processing_id", default_data_processing_id);
  put_opt(o, "default_source_file_id", default_source_file_id);
  return o;
}

/******************************************************************************/
RunMetadata::RunMetadata(const json::object& metadata)
    : raw_(metadata)
{
  if (const auto it = metadata.find("run");
      it != metadata.end() && it->value().is_object()) {
    run_ = Run::from_json(it->value().as_object());
  }

  if (const auto it = metadata.find("file_description");
      it != metadata.end() && it->value().is_object()) {
    file_description_ = FileDescription::from_json(it->value().as_object());
  }

  for (const auto& v : array_member(metadata, "software_list")) {
    if (v.is_object()) software_list_.push_back(Software::from_json(v.as_object()));
  }

  for (const auto& v : array_member(metadata, "instrument_configuration_list")) {
    if (v.is_object()) {
      instrument_configurations_.push_back(
          InstrumentConfiguration::from_json(v.as_object()));
    }
  }

  for (const auto& v : array_member(metadata, "data_processing_method_list")) {
    if (v.is_object()) {
      data_processings_.push_back(DataProcessing::from_json(v.as_object()));
    }
  }

  for (const auto& v : array_member(metadata, "sample_list")) {
    if (v.is_object()) samples_.push_back(Sample::from_json(v.as_object()));
  }
}

/******************************************************************************/
// WRT-2 — inverse of the parsing ctor: emit the typed blocks on top of the
// original raw_ object so unknown keys (e.g. scan_settings_list, cv_list)
// are preserved transparently; typed fields then override.
json::object RunMetadata::to_json() const
{
  // Start from raw_ so unknown blocks (e.g. scan_settings_list, cv_list)
  // are preserved transparently; typed fields then override.
  json::object o = raw_;

  if (run_.has_value()) o["run"] = run_->to_json();
  if (file_description_.has_value())
    o["file_description"] = file_description_->to_json();

  if (!software_list_.empty()) {
    json::array a;
    a.reserve(software_list_.size());
    for (const auto& s : software_list_)
      a.push_back(s.to_json());
    o["software_list"] = std::move(a);
  }

  if (!instrument_configurations_.empty()) {
    json::array a;
    a.reserve(instrument_configurations_.size());
    for (const auto& ic : instrument_configurations_)
      a.push_back(ic.to_json());
    o["instrument_configuration_list"] = std::move(a);
  }

  if (!data_processings_.empty()) {
    json::array a;
    a.reserve(data_processings_.size());
    for (const auto& dp : data_processings_)
      a.push_back(dp.to_json());
    o["data_processing_method_list"] = std::move(a);
  }

  if (!samples_.empty()) {
    json::array a;
    a.reserve(samples_.size());
    for (const auto& s : samples_)
      a.push_back(s.to_json());
    o["sample_list"] = std::move(a);
  }

  // cv_list is deliberately NOT injected here.  The index writer derives it
  // from every CV prefix actually present in the finished index -- including
  // whatever these run-metadata blocks carry -- so injecting an MS/UO-only list
  // here would shadow that derivation and hide a caller's other CVs.  A caller
  // that has already populated `cv_list` in the raw object is respected.
  return o;
}

} // namespace MzPeak
