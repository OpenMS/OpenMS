/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/schema/file.h"

#include "mzpeak/exception.h"

namespace MzPeak::Schema {

/******************************************************************************/
void parse_columns(const json::array& input, std::vector<File::Column>& output)
{
  auto get_string = [](const json::object& ob,
                       std::string_view key) -> std::optional<std::string> {
    auto it = ob.find(key);

    if (it != ob.end() && it->value().is_string()) {
      return std::string(it->value().as_string());
    } else {
      return {};
    }
  };

  auto get_bool = [](const json::object& ob,
                     std::string_view key) -> std::optional<bool> {
    auto it = ob.find(key);

    if (it != ob.end() && it->value().is_bool()) {
      return it->value().as_bool();
    } else {
      return {};
    }
  };

  output.reserve(input.size());

  for (const auto& column : input) {
    if (!column.is_object()) {
      // FIXME: emit a warning
      continue;
    }

    const auto& colobj(column.as_object());

    File::Column fc = {
        .name = get_string(colobj, "name").value_or(""),
        .path = get_string(colobj, "path")
                    .or_else(std::bind(get_string, colobj, "name"))
                    .value_or(""),
        .accession = get_string(colobj, "accession"),
        .unit = get_string(colobj, "unit"),
        .term_marker = get_bool(colobj, "term_marker").value_or(false),
    };

    // A mapping without a path binds nothing and would only ever produce a
    // lookup that silently misses; drop it at parse time instead.
    if (!fc.path.empty()) output.push_back(std::move(fc));
  }
}

/******************************************************************************/
File::File(const std::string& name)
    : file_name_(name)
{
}

/******************************************************************************/
File::File(const json::object& o)
try
    // Accept "path" (preferred) with "name" as backward-compatible fallback.
    : file_name_(o.contains("path") ? o.at("path").as_string()
                                    : o.at("name").as_string())
    , data_kind_(std::string_view(o.at("data_kind").as_string()))
    , entity_type_(std::string_view(o.at("entity_type").as_string()))
    , columns_() {
  auto cs = o.find("column_mapping");

  if (cs != o.end() && cs->value().is_array()) {
    parse_columns(cs->value().as_array(), columns_);
  }

  // Explicitly null is how an archive says "no digest", and is not an error.
  if (auto sum = o.find("checksum"); sum != o.end() && sum->value().is_string()) {
    checksum_ = std::string(sum->value().as_string());
  }
} catch (const std::exception& e) {
  throw MzPeak::JsonError(std::string("schema::File: ") + e.what());
}

/******************************************************************************/
std::optional<std::string> File::path_for(const std::string_view& accession) const
{
  for (const auto& c : columns_) {
    if (c.accession && *c.accession == accession) return c.path;
  }
  return std::nullopt;
}

/******************************************************************************/
std::optional<std::string> File::unit_for(const std::string_view& accession) const
{
  for (const auto& c : columns_) {
    if (c.accession && *c.accession == accession) return c.unit;
  }
  return std::nullopt;
}

/******************************************************************************/
bool File::operator==(const File& other) const
{
  return file_name_ == other.file_name_;
};
} // namespace MzPeak::Schema
