/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <boost/json.hpp>
#include <optional>
#include <string>
#include <vector>

namespace MzPeak {

/**
 * RDR-24 — typed file-level (run-level) metadata blocks parsed from the
 * `metadata{}` object of `mzpeak_index.json`.
 *
 * The blocks are heterogeneous CV-param trees.  This model gives the
 * highest-value fields (block ids, names, references) named accessors while
 * retaining every CV term as a flat `CvParam`.  It is deliberately *partial*:
 * a `CvParam` keeps `accession`/`name`/`unit`/`value` but does not interpret
 * the term or resolve its CV; callers that need more reach for the raw
 * Boost.JSON value preserved alongside the typed blocks
 * (`RunMetadata::raw()`).
 *
 * WRT-2 — each block and `RunMetadata` also provides a `to_json()` serializer,
 * the inverse of the `from_json` parsers, so the C++ writer can emit these
 * run-level blocks into the index `metadata{}` object.  The serialization
 * mirrors the parser's field names, so `parse(serialize(x)) == x` holds for
 * the typed fields.
 */

/**
 * A single controlled-vocabulary parameter, mirroring the
 * `{accession, name, unit, value}` shape used throughout the index metadata.
 * Any of the fields may be absent (JSON `null`); `value` is captured as its
 * stringified form (the source mixes strings and numbers).
 */
struct CvParam {
  std::optional<std::string> accession;
  std::optional<std::string> name;
  std::optional<std::string> unit;
  std::optional<std::string> value;

  /// Parse a CvParam from its JSON parameter object.
  static CvParam from_json(const boost::json::object&);

  /// Serialize to its `{accession, name, unit, value}` JSON object (the
  /// inverse of from_json); absent optionals are omitted.
  boost::json::object to_json() const;
};

/// Parse a JSON array of parameter objects into a `CvParam` vector.
std::vector<CvParam> cv_params_from_json(const boost::json::array&);

/// Serialize a `CvParam` vector to a JSON array of parameter objects.
boost::json::array cv_params_to_json(const std::vector<CvParam>&);

/**
 * A software package (`software_list[]` entry): id, version, and CV params
 * (e.g. the software-name term).
 */
struct Software {
  std::string id;
  std::optional<std::string> version;
  std::vector<CvParam> parameters;

  /// Parse a Software from its JSON object representation.
  static Software from_json(const boost::json::object&);

  /// Serialize to its JSON object representation (inverse of from_json).
  boost::json::object to_json() const;
};

/**
 * One ordered component of an instrument configuration (source / analyzer /
 * detector), with its CV params.
 */
struct Component {
  std::optional<std::string> component_type;
  std::optional<std::int64_t> order;
  std::vector<CvParam> parameters;

  /// Parse a Component from its JSON object representation.
  static Component from_json(const boost::json::object&);

  /// Serialize to its JSON object representation (inverse of from_json).
  boost::json::object to_json() const;
};

/**
 * An instrument configuration (`instrument_configuration_list[]` entry):
 * numeric id, top-level CV params (instrument model, serial number, …), an
 * optional software reference, and the ordered component list.
 */
struct InstrumentConfiguration {
  std::optional<std::int64_t> id;
  std::optional<std::string> software_reference;
  std::vector<CvParam> parameters;
  std::vector<Component> components;

  /// Parse a InstrumentConfiguration from its JSON object representation.
  static InstrumentConfiguration from_json(const boost::json::object&);

  /// Serialize to its JSON object representation (inverse of from_json).
  boost::json::object to_json() const;
};

/**
 * One ordered method within a data-processing chain, with the software that
 * performed it and its CV params.
 */
struct ProcessingMethod {
  std::optional<std::int64_t> order;
  std::optional<std::string> software_reference;
  std::vector<CvParam> parameters;

  /// Parse a ProcessingMethod from its JSON object representation.
  static ProcessingMethod from_json(const boost::json::object&);

  /// Serialize to its JSON object representation (inverse of from_json).
  boost::json::object to_json() const;
};

/**
 * A data-processing chain (`data_processing_method_list[]` entry): a string id
 * and its ordered methods.
 */
struct DataProcessing {
  std::string id;
  std::vector<ProcessingMethod> methods;

  /// Parse a DataProcessing from its JSON object representation.
  static DataProcessing from_json(const boost::json::object&);

  /// Serialize to its JSON object representation (inverse of from_json).
  boost::json::object to_json() const;
};

/**
 * A sample (`sample_list[]` entry): id, name, and CV params.
 */
struct Sample {
  std::string id;
  std::optional<std::string> name;
  std::vector<CvParam> parameters;

  /// Parse a Sample from its JSON object representation.
  static Sample from_json(const boost::json::object&);

  /// Serialize to its JSON object representation (inverse of from_json).
  boost::json::object to_json() const;
};

/**
 * A source file (`file_description.source_files[]` entry): id, name, location,
 * and CV params (checksum, native-id format, …).
 */
struct SourceFile {
  std::string id;
  std::optional<std::string> name;
  std::optional<std::string> location;
  std::vector<CvParam> parameters;

  /// Parse a SourceFile from its JSON object representation.
  static SourceFile from_json(const boost::json::object&);

  /// Serialize to its JSON object representation (inverse of from_json).
  boost::json::object to_json() const;
};

/**
 * The file description (`file_description` block): the run's content CV params
 * (spectrum-type terms) and the list of source files.
 */
struct FileDescription {
  std::vector<CvParam> contents;
  std::vector<SourceFile> source_files;

  /// Parse a FileDescription from its JSON object representation.
  static FileDescription from_json(const boost::json::object&);

  /// Serialize to its JSON object representation (inverse of from_json).
  boost::json::object to_json() const;
};

/**
 * The run block (`run`): run id, acquisition start timestamp, and the default
 * instrument / data-processing / source-file references that downstream
 * OpenMS `ExperimentalSettings` mapping (RDR-19) hangs off of.
 */
struct Run {
  std::optional<std::string> id;
  std::optional<std::string> start_time;
  std::optional<std::int64_t> default_instrument_id;
  std::optional<std::string> default_data_processing_id;
  std::optional<std::string> default_source_file_id;

  /// Parse a Run from its JSON object representation.
  static Run from_json(const boost::json::object&);

  /// Serialize to its JSON object representation (inverse of from_json).
  boost::json::object to_json() const;
};

/**
 * The complete parsed run-level metadata.  Each block is typed; the original
 * `metadata{}` object is retained verbatim (`raw()`) for callers that need a
 * term this partial model does not surface.  `scan_settings_list` is kept only
 * in raw form (it is empty in every reference file and has no fixed shape).
 */
class RunMetadata {
public:
  RunMetadata() = default;

  /// Parse the index `metadata{}` object (everything except `version`).
  explicit RunMetadata(const boost::json::object&);

  /// The run block, if a `run` object was present.
  const std::optional<Run>& run() const { return run_; }

  /// The file description, if a `file_description` object was present.
  const std::optional<FileDescription>& file_description() const
  {
    return file_description_;
  }

  /// `software_list` (empty if absent).
  const std::vector<Software>& software_list() const { return software_list_; }

  /// `instrument_configuration_list` (empty if absent).
  const std::vector<InstrumentConfiguration>& instrument_configurations() const
  {
    return instrument_configurations_;
  }

  /// `data_processing_method_list` (empty if absent).
  const std::vector<DataProcessing>& data_processings() const
  {
    return data_processings_;
  }

  /// `sample_list` (empty if absent).
  const std::vector<Sample>& samples() const { return samples_; }

  /// `true` if the index carried no run-level metadata block.  Tests the typed
  /// blocks rather than `raw_`, so a `metadata{}` that holds only the format
  /// `version` (handled separately by Index) still reports empty.
  bool empty() const
  {
    return !run_.has_value() && !file_description_.has_value() &&
           software_list_.empty() && instrument_configurations_.empty() &&
           data_processings_.empty() && samples_.empty();
  }

  /// The verbatim `metadata{}` object, for blocks/terms not typed above.
  const boost::json::object& raw() const { return raw_; }

  /// WRT-2 — serialize the typed blocks into a `metadata{}` JSON object (the
  /// inverse of the parsing ctor), suitable for merging into the index's
  /// emitted `metadata{}` alongside `version`.  Only the typed blocks are
  /// emitted; absent optionals and empty lists are omitted, so
  /// `RunMetadata(x.to_json())` reproduces the typed fields of `x`.
  boost::json::object to_json() const;

private:
  boost::json::object raw_;
  std::optional<Run> run_;
  std::optional<FileDescription> file_description_;
  std::vector<Software> software_list_;
  std::vector<InstrumentConfiguration> instrument_configurations_;
  std::vector<DataProcessing> data_processings_;
  std::vector<Sample> samples_;
};

} // namespace MzPeak
