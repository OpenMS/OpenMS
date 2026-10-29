/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/util/json_writer.h"

#include <algorithm>
#include <boost/json.hpp>
#include <cctype>
#include <cstdio>
#include <set>
#include <string_view>

namespace MzPeak::Util {

namespace json = boost::json;

/******************************************************************************/
std::string point_spectra_array_index_json()
{
  // Mirrors the `spectrum_array_index` blob stored in the reference
  // file: the index column (point.spectrum_index) is synthesized by the
  // reader and so only the m/z and intensity entries are emitted here.
  json::object root;
  root["prefix"] = "point";

  json::array entries;

  json::object mz;
  mz["context"] = "spectrum";
  mz["path"] = "point.mz";
  mz["data_type"] = "MS:1000523";  // 64-bit float
  mz["array_type"] = "MS:1000514"; // m/z array
  mz["array_name"] = "m/z array";
  mz["unit"] = "MS:1000040"; // m/z
  mz["buffer_format"] = "point";
  mz["transform"] = nullptr;
  mz["data_processing_id"] = nullptr;
  mz["buffer_priority"] = "primary";
  mz["sorting_rank"] = 0;
  entries.push_back(std::move(mz));

  json::object intensity;
  intensity["context"] = "spectrum";
  intensity["path"] = "point.intensity";
  intensity["data_type"] = "MS:1000521";  // 32-bit float
  intensity["array_type"] = "MS:1000515"; // intensity array
  intensity["array_name"] = "intensity array";
  intensity["unit"] = "MS:1000131"; // number of detector counts
  intensity["buffer_format"] = "point";
  intensity["transform"] = nullptr;
  intensity["data_processing_id"] = nullptr;
  intensity["buffer_priority"] = "primary";
  intensity["sorting_rank"] = nullptr;
  entries.push_back(std::move(intensity));

  root["entries"] = std::move(entries);

  return json::serialize(root);
}

/******************************************************************************/
namespace {

/// The controlled vocabularies whose CURIEs this writer emits.
///
/// conformance.md requires an archive to declare every CV prefix it uses, with
/// a uri and a version.  Nothing enforces it and the reference fixtures omit
/// the list entirely, but a consumer that cannot resolve "MS:1000235" has no
/// way to learn what the term means.
/// A controlled vocabulary this writer can pin to a version and a URI.
struct KnownCv {
  const char* id;
  const char* full_name;
  const char* version;
  const char* uri;
};

/// The vocabularies a mass-spectrometry archive plausibly cites.  MS and UO are
/// always present (this writer's array index and column mappings use them); the
/// rest cover terms a caller may carry through run metadata.
const KnownCv kKnownCvs[] = {
    {"MS", "PSI-MS controlled vocabulary", "4.1.209",
     "https://raw.githubusercontent.com/HUPO-PSI/psi-ms-CV/master/psi-ms.obo"},
    {"UO", "Units of Measurement Ontology", "releases/2023-05-25",
     "http://purl.obolibrary.org/obo/uo.owl"},
    {"PATO", "Phenotype And Trait Ontology", "releases/2023-05-18",
     "http://purl.obolibrary.org/obo/pato.owl"},
    {"NCIT", "NCI Thesaurus", "23.03e", "http://purl.obolibrary.org/obo/ncit.owl"},
    {"UNIMOD", "Unimod protein modifications", "2023-05",
     "http://www.unimod.org/obo/unimod.obo"},
    {"MOD", "PSI-MOD protein modification ontology", "1.031.4",
     "http://purl.obolibrary.org/obo/mod.obo"},
    {"BTO", "BRENDA Tissue Ontology", "2021-10-26",
     "http://purl.obolibrary.org/obo/bto.owl"},
    {"CL", "Cell Ontology", "2023-04-20", "http://purl.obolibrary.org/obo/cl.owl"},
    {"NCBITaxon", "NCBI organismal taxonomy", "2023-06-20",
     "http://purl.obolibrary.org/obo/ncbitaxon.owl"},
};

/// Collect the CURIE prefixes ("MS" from "MS:1000235") of every accession-shaped
/// string reachable from @p value.  A conformant archive must declare every one.
void collect_cv_prefixes(const json::value& value, std::set<std::string>& out)
{
  if (value.is_object()) {
    for (const auto& [key, v] : value.as_object())
      collect_cv_prefixes(v, out);
  } else if (value.is_array()) {
    for (const auto& v : value.as_array())
      collect_cv_prefixes(v, out);
  } else if (value.is_string()) {
    std::string_view s(value.as_string());
    // A CURIE is "<prefix>:<local>" where the prefix is a CV id.  Match
    // conservatively: an all-alphanumeric prefix, a colon, then a non-empty
    // local part, and no whitespace or second colon.
    auto colon = s.find(':');
    if (colon == std::string_view::npos || colon == 0 || colon + 1 >= s.size()) {
      return;
    }
    std::string_view prefix = s.substr(0, colon);
    std::string_view local = s.substr(colon + 1);
    auto ok = [](std::string_view p) {
      return !p.empty() &&
             std::ranges::all_of(p, [](unsigned char c) { return std::isalnum(c); });
    };
    if (ok(prefix) && ok(local)) out.emplace(prefix);
  }
}

/// Build `cv_list` from the vocabularies the archive actually uses.
///
/// The previous version hardcoded MS and UO regardless of content, so an
/// archive carrying, say, an NCIT term in its run metadata declared a CV list
/// that did not mention NCIT -- a conformance failure the spec calls out
/// (`cv_list` MUST name every prefix used).  MS and UO are still always present
/// because the array index and column mappings use them; everything else is
/// derived from what is actually written.
json::array cv_list_for(const json::object& index_root)
{
  std::set<std::string> prefixes{"MS", "UO"};
  collect_cv_prefixes(index_root, prefixes);

  json::array list;
  for (const auto& prefix : prefixes) {
    const KnownCv* known = nullptr;
    for (const auto& cv : kKnownCvs) {
      if (prefix == cv.id) {
        known = &cv;
        break;
      }
    }

    json::object entry;
    entry["id"] = prefix;
    if (known != nullptr) {
      entry["full_name"] = known->full_name;
      entry["version"] = known->version;
      entry["uri"] = known->uri;
    } else {
      // A prefix this writer cannot pin.  Declaring it -- even without an
      // authoritative version -- keeps every used prefix named, which is what
      // the spec requires; the placeholder is visible rather than silent.
      std::fprintf(stderr,
                   "mzpeak: cv_list: unknown CV prefix '%s'; declaring it with a "
                   "placeholder version/uri\n",
                   prefix.c_str());
      entry["full_name"] = prefix;
      entry["version"] = "unknown";
      entry["uri"] = "urn:mzpeak:unknown-cv:" + prefix;
    }
    list.push_back(std::move(entry));
  }
  return list;
}

// Build the shared `files` array + `metadata{}` skeleton.  `metadata` is
// supplied pre-populated with the run-level blocks (empty for the no-metadata
// overload); `version` is then stamped on, always overriding any pre-existing
// `version` member.
std::string index_json(const std::vector<IndexFileEntry>& files,
                       const std::string& version,
                       json::object metadata)
{
  json::object root;

  json::array file_array;
  file_array.reserve(files.size());

  for (const auto& file : files) {
    json::object o;
    o["name"] = file.name;
    o["entity_type"] = file.entity_type;
    o["data_kind"] = file.data_kind;
    // The reference reader resolves metadata columns through this mapping; a
    // column missing from it is reported as "unspecified column" and ignored.
    json::array mapping;
    for (const auto& m : file.column_mapping) {
      json::object e;
      e["name"] = m.name;
      e["path"] = m.path;
      e["accession"] = m.accession;
      if (m.unit.empty()) {
        e["unit"] = nullptr;
      } else {
        e["unit"] = m.unit;
      }
      mapping.push_back(std::move(e));
    }
    o["column_mapping"] = std::move(mapping);
    o["parameters"] = json::array();
    if (file.checksum.empty()) {
      o["checksum"] = nullptr;
    } else {
      o["checksum"] = file.checksum;
    }
    file_array.push_back(std::move(o));
  }

  root["files"] = std::move(file_array);

  metadata["version"] = version;
  root["metadata"] = std::move(metadata);

  // Derive cv_list from every CURIE prefix reachable in the finished index --
  // files, column mappings and run metadata alike -- so the archive names every
  // CV it uses.  A caller-supplied cv_list is respected and not overwritten.
  if (!root["metadata"].as_object().contains("cv_list")) {
    root["metadata"].as_object()["cv_list"] = cv_list_for(root);
  }

  return json::serialize(root);
}

} // namespace

/******************************************************************************/
std::string mzpeak_index_json(const std::vector<IndexFileEntry>& files,
                              const std::string& version)
{
  return index_json(files, version, json::object{});
}

/******************************************************************************/
// WRT-2 — merge the serialized run-level metadata blocks into `metadata{}`.
std::string mzpeak_index_json(const std::vector<IndexFileEntry>& files,
                              const std::string& version,
                              const json::object& run_metadata)
{
  return index_json(files, version, run_metadata);
}

/******************************************************************************/
namespace {

/// One array-index entry.  Every field the schema requires is written, including
/// `unit`: the reference's own serializable type models it as optional, but the
/// JSON Schema requires a CURIE and a canonical array always has a real unit.
json::object array_entry(const char* context,
                         const char* path,
                         const char* data_type,
                         const char* array_type,
                         const char* array_name,
                         const std::string& unit,
                         bool primary,
                         bool sorted)
{
  json::object e;
  e["context"] = context;
  e["path"] = path;
  e["data_type"] = data_type;
  e["array_type"] = array_type;
  e["array_name"] = array_name;
  e["unit"] = unit;
  e["buffer_format"] = "point";
  e["transform"] = nullptr;
  e["data_processing_id"] = nullptr;
  e["buffer_priority"] = primary ? json::value("primary") : json::value(nullptr);
  e["sorting_rank"] = sorted ? json::value(0) : json::value(nullptr);
  return e;
}

} // namespace

/******************************************************************************/
std::string point_chromatograms_array_index_json(const std::string& intensity_unit)
{
  json::object root;
  root["prefix"] = "point";

  json::array entries;
  // Time is stored in MINUTES (UO:0000031), which is what the specification
  // recommends and what every reference file uses.  The reader converts to
  // seconds using this declaration rather than assuming.
  entries.push_back(array_entry("chromatogram", "point.time", "MS:1000523",
                                "MS:1000595", "time array", "UO:0000031",
                                /*primary=*/true, /*sorted=*/true));
  entries.push_back(array_entry("chromatogram", "point.intensity", "MS:1000521",
                                "MS:1000515", "intensity array", intensity_unit,
                                /*primary=*/true, /*sorted=*/false));
  root["entries"] = std::move(entries);
  return json::serialize(root);
}

/******************************************************************************/
std::string point_wavelength_array_index_json(const std::string& intensity_unit)
{
  json::object root;
  root["prefix"] = "point";

  json::array entries;
  entries.push_back(array_entry("wavelength_spectrum", "point.wavelength",
                                "MS:1000521", "MS:1000617", "wavelength array",
                                "UO:0000018", /*primary=*/true, /*sorted=*/true));
  entries.push_back(array_entry("wavelength_spectrum", "point.intensity",
                                "MS:1000521", "MS:1000515", "intensity array",
                                intensity_unit, /*primary=*/true,
                                /*sorted=*/false));
  root["entries"] = std::move(entries);
  return json::serialize(root);
}

} // namespace MzPeak::Util
