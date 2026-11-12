/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/metadata/table.h"

#include <parquet/api/reader.h>

#include "mzpeak/exception.h"
#include "mzpeak/schema/group.h"
#include "mzpeak/util/metadata_model.h"
#include "mzpeak/util/parquet.h"

namespace MzPeak::Metadata {

/******************************************************************************/
struct Table::Impl {
  Impl(std::unique_ptr<Util::Parquet> parquet);
  ~Impl();

  std::unique_ptr<Util::Parquet> parquet_;
  std::unique_ptr<Util::Parquet> scans_;
  std::unique_ptr<Util::Parquet> precursors_;
  std::unique_ptr<Util::Parquet> selected_ions_;
  std::optional<std::size_t> n_entries;
  std::string index_field_name_;
};

/******************************************************************************/
Table::Impl::Impl(std::unique_ptr<Util::Parquet> parquet)
    : parquet_(std::move(parquet))
    , index_field_name_("index")
{
  auto file = parquet_->index_file();

  if (!file.data_kind().is_metadata()) {
    std::string msg("file is not a metadata file: " + file.file_name());
    throw InvalidFormatError(msg);
  }

  auto type = file.data_kind().type();
  if (type.has_value()) {
    using enum Schema::DataKind::Type;

    switch (type.value()) {
    case DataArray:
    case Peaks:
    case Metadata:
    case Proprietary:
      break;
    case Scans:
    case Precursors:
    case SelectedIons:
    case Products:
      index_field_name_ = "source_index";
      break;
    }
  }
}

/******************************************************************************/
Table::Impl::~Impl() = default;

/******************************************************************************/
Table::Table(std::unique_ptr<Util::Parquet> parquet)
    : impl_(std::make_unique<Impl>(std::move(parquet)))
{
}

/******************************************************************************/
Table::~Table() = default;

/******************************************************************************/
std::shared_ptr<Schema::Group> Table::group(std::string_view name) const
{
  const std::shared_ptr<Schema::GroupMap>& map = impl_->parquet_->groups();
  auto it = map->find(std::string(name));

  if (it == map->end()) {
    throw InvalidFormatError("schema is missing the " + std::string(name) +
                             " group");
  } else {
    return it->second;
  }
}

/******************************************************************************/
std::unique_ptr<Util::Slice>
Table::indexed(uint64_t index,
               const std::shared_ptr<Schema::Group>& group,
               const Util::Projection& projection) const
{
  auto index_field = group->field(impl_->index_field_name_);

  if (!index_field.has_value()) {
    throw InvalidFormatError("metadata column missing: " + impl_->index_field_name_);
  }

  auto column = std::make_pair(group, index_field.value());
  Util::Query q = Util::Query::Builder(column).eq(index);

  auto plan = impl_->parquet_->planner(q).plan();
  return impl_->parquet_->executor(projection).execute(plan);
}

/******************************************************************************/
void Table::add_facet(Schema::DataKind::Type kind, std::unique_ptr<Util::Parquet> p)
{
  using enum Schema::DataKind::Type;
  switch (kind) {
  case Scans:
    impl_->scans_ = std::move(p);
    break;
  case Precursors:
    impl_->precursors_ = std::move(p);
    break;
  case SelectedIons:
    impl_->selected_ions_ = std::move(p);
    break;
  default:
    break; // not a facet; ignore
  }
}

/******************************************************************************/
Util::IndexMap<SpectrumMetadata>
Table::read_spectrum_metadata(MetadataDetail detail) const
{
  Util::SpectraMetadataFiles files;
  files.primary = impl_->parquet_.get();
  files.scans = impl_->scans_.get();
  files.precursors = impl_->precursors_.get();
  files.selected_ions = impl_->selected_ions_.get();
  return Util::read_spectra_metadata(files, detail);
}

/******************************************************************************/
std::optional<Util::IndexMap<MinimalSpectrumMetadata>>
Table::read_minimal_spectrum_metadata() const
{
  Util::SpectraMetadataFiles files;
  files.primary = impl_->parquet_.get();
  files.scans = impl_->scans_.get();
  files.precursors = impl_->precursors_.get();
  files.selected_ions = impl_->selected_ions_.get();
  return Util::read_minimal_spectra_metadata(files);
}

} // namespace MzPeak::Metadata
