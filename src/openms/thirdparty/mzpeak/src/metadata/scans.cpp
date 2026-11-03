/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/metadata/scans.h"
#include "mzpeak/metadata/table.h"
#include "mzpeak/schema/group.h"
#include "mzpeak/util/projection.h"
#include "mzpeak/util/slice.h"

namespace MzPeak::Metadata {

using namespace MzPeak::Schema;
using namespace MzPeak::Util;

/******************************************************************************/
Scans::Scans()
    : raw_()
{
}

/******************************************************************************/
Scans::Scans(std::unique_ptr<Util::Parquet> parquet, uint64_t index)
    : raw_()
{
  Table table(std::move(parquet));
  std::shared_ptr<Schema::Group> group = table.group("root");

  Projection projection;
  auto scan_time_field = projection.project(group, Group::CVType("MS", "1000016"));

  std::unique_ptr<Slice> slice = table.indexed(index, group, projection);
  slice->non_null(scan_time_field, raw_.scan_start_time);
}

} // namespace MzPeak::Metadata
