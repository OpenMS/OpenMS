/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/util/manager.h"

#include "mzpeak/exception.h"

namespace MzPeak::Util {

/******************************************************************************/
const static char* INDEX_FILE_NAME = "mzpeak_index.json";

/******************************************************************************/
namespace json = boost::json;

/******************************************************************************/
void parse_index(std::shared_ptr<MzPeak::IO::Archive>& archive,
                 std::vector<Schema::File>& files,
                 std::string& version,
                 ImsCalibration& ims,
                 RunMetadata& metadata)
{
  auto file = archive->read_file(INDEX_FILE_NAME);
  uint8_t buffer[64 * 1024];
  std::optional<std::size_t> bytes;

  json::stream_parser parser;
  boost::system::error_code ec;

  do {
    bytes = file->read(buffer, sizeof(buffer));

    if (bytes.has_value() && *bytes > 0) {
      parser.write(reinterpret_cast<char const*>(buffer), *bytes, ec);
    }
  } while (bytes.has_value() && !ec);

  if (!ec) parser.finish(ec);
  if (ec) throw MzPeak::JsonError(ec.message());

  json::value v = parser.release();
  json::object o = v.as_object();

  if (const auto it = o.find("files"); it != o.end() && it->value().is_array()) {
    const json::array file_list(it->value().as_array());
    files.reserve(file_list.size());

    for (const auto& file_obj : file_list) {
      if (file_obj.is_object()) {
        files.push_back(Schema::File(file_obj.as_object()));
      }
    }
  }

  // Format version from metadata.version (e.g. "0.9.0").
  if (const auto it = o.find("metadata"); it != o.end() && it->value().is_object()) {
    const json::object& meta = it->value().as_object();

    // The typed run-level blocks (run, file_description, software_list, ...).
    metadata = RunMetadata(meta);
    if (const auto v = meta.find("version");
        v != meta.end() && v->value().is_string()) {
      version = v->value().as_string().c_str();
    }

    // ims-compact: {"a": ..., "b": ..., "mz_from_tof": "(a + b*tof)^2"}.
    // Writers emit this either as a JSON object or as a STRING holding that
    // object; accept both rather than silently reading no calibration, which
    // would make every spectrum of such a file come back empty.
    if (const auto ic = meta.find("ims_calibration"); ic != meta.end()) {
      json::value parsed;
      if (ic->value().is_string()) {
        boost::system::error_code pe;
        parsed = json::parse(ic->value().as_string(), pe);
        if (pe) parsed = nullptr;
      }
      const json::value& calv = ic->value().is_string() ? parsed : ic->value();

      if (calv.is_object()) {
        const json::object& cal = calv.as_object();
        const auto number = [&cal](const char* key, double& dst) {
          const auto f = cal.find(key);
          if (f == cal.end() || !f->value().is_number()) return false;
          dst = f->value().to_number<double>();
          return true;
        };
        // The archive states which transform its coefficients belong to.  We
        // implement exactly one, so a file naming a different one must be
        // refused rather than run through this one: the coefficients would be
        // consumed by the wrong formula and every m/z would be wrong while
        // remaining entirely plausible.
        //
        // Verified against a real Bruker timsTOF archive, whose calibration
        // reads {"a": 9.9995…, "b": 4.9255e-05, "mz_from_tof": "(a + b*tof)^2",
        // "tof_encoding": "absolute", "codec": "ims-compact"}.
        const auto declares = [&cal](const char* key, std::string_view expected) {
          const auto f = cal.find(key);
          if (f == cal.end() || !f->value().is_string()) return true; // not stated
          return std::string_view(f->value().as_string()) == expected;
        };

        // Both coefficients are required: half a calibration is not one.
        const bool have_a = number("a", ims.a);
        const bool have_b = number("b", ims.b);
        ims.valid = have_a && have_b;

        // The declaration guards below apply only when the block actually
        // hands us coefficients to apply.  A block with no `a`/`b` is
        // DESCRIPTIVE -- the grid-transform layout writes one to document its
        // models and column names, and there is no formula there for this
        // reader to get wrong.  Refusing it would reject a well-formed archive
        // over metadata nothing here consumes, which is what happened before
        // this ordering was fixed.
        if (ims.valid) {
          if (!declares("mz_from_tof", "(a + b*tof)^2")) {
            throw MzPeak::JsonError(
                "ims_calibration declares an m/z transform this reader does not "
                "implement; refusing rather than applying the wrong one");
          }
          if (!declares("tof_encoding", "absolute")) {
            throw MzPeak::JsonError(
                "ims_calibration declares a non-absolute TOF encoding; this "
                "reader would treat the stored values as absolute");
          }
        }
      }
    }
  }

  // Reject a future, potentially incompatible MAJOR version (the spec is a
  // pre-1.0 living standard; all 0.x are accepted).
  if (!version.empty()) {
    std::string major(version.substr(0, version.find('.')));
    if (major != "0") {
      throw MzPeak::JsonError("unsupported mzPeak format version '" + version +
                              "' (this reader supports 0.x)");
    }
  }
}

/******************************************************************************/
Manager::Manager(std::unique_ptr<MzPeak::IO::Archive> archive)
    : archive_(std::move(archive))
    , files_()
    , version_()
    , ims_()
    , metadata_()
{
  parse_index(archive_, files_, version_, ims_, metadata_);
}

/******************************************************************************/
const std::vector<Schema::File>& Manager::files() const { return files_; }

/******************************************************************************/
std::vector<Schema::File>::const_iterator
Manager::find_file(std::string_view name) const
{
  return std::ranges::find(files_, name, &Schema::File::file_name);
}

/******************************************************************************/
std::vector<Schema::File>::const_iterator
Manager::find_file(Schema::EntityType::Type et, Schema::DataKind::Type dkt) const
{
  return std::ranges::find_if(files_, [&et, &dkt](const auto& file) -> bool {
    auto et_type = file.entity_type().type();
    auto dk_type = file.data_kind().type();

    return et_type.has_value() && et_type.value() == et && dk_type.has_value() &&
           dk_type.value() == dkt;
  });
}

/******************************************************************************/
std::shared_ptr<const Manager::SpectrumMetadataMap>
Manager::spectrum_metadata(MetadataDetail detail,
                           const std::function<SpectrumMetadataMap()>& build) const
{
  // The lock is held across the BUILD, not just the lookup.  Two threads
  // opening spectra at the same moment would otherwise both read the whole
  // metadata table -- the exact duplication this cache exists to remove -- and
  // the loser's copy would then be thrown away.  Waiting is cheaper than
  // building, and the second thread finds it cached.
  std::lock_guard<std::mutex> guard(spectrum_metadata_mutex_);

  auto it = spectrum_metadata_.find(detail);
  if (it != spectrum_metadata_.end()) return it->second;

  auto map = std::make_shared<const SpectrumMetadataMap>(build());
  spectrum_metadata_.emplace(detail, map);
  return map;
}

/******************************************************************************/
std::shared_ptr<const Manager::MinimalSpectrumMetadataMap>
Manager::minimal_spectrum_metadata(
    const std::function<std::optional<MinimalSpectrumMetadataMap>()>& build) const
{
  // Locked across the build for the reason spectrum_metadata() gives.
  std::lock_guard<std::mutex> guard(spectrum_metadata_mutex_);
  if (!minimal_spectrum_metadata_) {
    auto map = build();
    minimal_spectrum_metadata_ =
        map ? std::make_shared<const MinimalSpectrumMetadataMap>(std::move(*map))
            : nullptr;
  }
  return *minimal_spectrum_metadata_;
}

/******************************************************************************/
std::unique_ptr<Util::Parquet> Manager::parquet(const Schema::File& file) const
{
  // Every Parquet opens its OWN handle -- the file position is per reader and
  // cannot be shared -- but the footer is the same for every handle on one
  // member, and reading and parsing it was most of opening one.  A per-thread
  // Spectra over a 717,924-spectrum run took ~2.5 ms to open its six members,
  // 1.0 ms of it for the 341-row-group peaks footer alone, and kept a parsed
  // copy of every footer.  So the first Parquet over a member parses it and
  // every later one reuses that: ~1 ms per Spectra now.
  //
  // Sharing is safe because a parsed FileMetaData is only ever read: every
  // accessor builds its row-group and column-chunk views afresh from it, and
  // Arrow's own dataset scanner hands one to concurrent readers the same way.
  // The lock covers only this map.  Two threads that open the same member for
  // the first time at once both parse it; the first result is kept, and the
  // other lives only as long as its reader.
  //
  // Only the footer.  The StatsIndex built from it stays per Parquet: its
  // sorted-column lookups take its mutex on every query and its group hint is
  // per-reader state, so sharing it would turn an uncontended lock per thread
  // into one that every reader thread contends for.
  Parquet::file_metadata_t footer;
  {
    std::lock_guard<std::mutex> guard(footers_mutex_);
    if (auto it = footers_.find(file.file_name()); it != footers_.end())
      footer = it->second;
  }

  std::unique_ptr<IO::File> data(archive_->read_file(file.file_name()));
  auto parquet = std::make_unique<Util::Parquet>(std::move(data), file,
                                                 row_group_cache_, footer);
  if (!footer) {
    std::lock_guard<std::mutex> guard(footers_mutex_);
    footers_.emplace(file.file_name(), parquet->file_metadata());
  }
  return parquet;
}

/******************************************************************************/
std::unique_ptr<MzPeak::IO::File> Manager::read_member(const std::string& name) const
{
  return archive_->read_file(name);
}

} // namespace MzPeak::Util
