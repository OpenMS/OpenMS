/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#include "mzpeak/schema/group.h"
#include "mzpeak/util/projection.h"

namespace MzPeak::Util {

/******************************************************************************/
Projection::Projection(std::size_t reserve)
    : projections_({})
{
  projections_.reserve(reserve);
}

/******************************************************************************/
Projection::Result Projection::project(const Schema::Column& c)
{
  projections_.push_back(c);
  return c;
}

/******************************************************************************/
Projection::Result Projection::project(
    const std::shared_ptr<Schema::Group>& group,
    const std::optional<std::shared_ptr<const Schema::Group::Field>>& field)
{
  if (field.has_value()) {
    auto col = std::make_pair(group, field.value());
    projections_.push_back(col);
    return col;
  } else {
    return {};
  }
}

/******************************************************************************/
Projection::Result Projection::project(const std::shared_ptr<Schema::Group>& group,
                                       std::string_view name)
{
  return project(group, group->field(std::move(name)));
}

/******************************************************************************/
Projection::Result Projection::project(const std::shared_ptr<Schema::Group>& group,
                                       Schema::Group::CVType&& cvt)
{
  return project(group, group->field(std::move(cvt)));
}

/******************************************************************************/
const std::vector<Schema::Column>& Projection::get() const { return projections_; }

} // namespace MzPeak::Util
