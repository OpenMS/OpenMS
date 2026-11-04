/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include <memory>

#include "mzpeak/metadata/spectrum.h"
#include "mzpeak/metadata/table.h"
#include "mzpeak/schema/group.h"
#include "mzpeak/util/projection.h"
#include "mzpeak/util/slice.h"

namespace MzPeak::Metadata {

using namespace MzPeak::Schema;
using namespace MzPeak::Util;

/******************************************************************************/
Spectrum::Spectrum(std::unique_ptr<Util::Parquet> parquet, uint64_t index)
    : ms_level_()
    , scan_time_()
    , delta_model_()
{
  Table table(std::move(parquet));
  std::shared_ptr<Schema::Group> group = table.group("root");

  Projection projection;
  auto level_field = projection.project(group, Group::CVType("MS", "1000511"));
  auto scan_time_field = projection.project(group, "time");
  auto delta_field = projection.project(group, "mz_delta_model");

  std::unique_ptr<Slice> slice = table.indexed(index, group, projection);

  slice->scalar(level_field, ms_level_);
  slice->scalar(scan_time_field, scan_time_);
  slice->list(delta_field, delta_model_);
}

} // namespace MzPeak::Metadata
