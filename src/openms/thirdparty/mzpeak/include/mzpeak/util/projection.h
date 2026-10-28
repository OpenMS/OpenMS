/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <string_view>
#include <vector>

#include "mzpeak/schema/group.h"

namespace MzPeak::Util {

/**
 * Columns that should be projected out of an executed query.
 */
class Projection final {
public:
  /// Result type.
  using Result = std::optional<Schema::Column>;

  /// Constructor.
  explicit Projection(std::size_t reserve = 0);

  /// Destructor.
  ~Projection() = default;

  /**
   * Project a column directly.
   */
  Result project(const Schema::Column&);

  /**
   * Project a column if the field exists.
   */
  Result project(const std::shared_ptr<Schema::Group>&,
                 const std::optional<std::shared_ptr<const Schema::Group::Field>>&);

  /**
   * Project a column using a field name.
   */
  Result project(const std::shared_ptr<Schema::Group>&, std::string_view);

  /**
   * Look up a CV type and project that.
   */
  Result project(const std::shared_ptr<Schema::Group>&, Schema::Group::CVType&&);

  /**
   * Return the vector of projected columns.
   */
  const std::vector<Schema::Column>& get() const;

private:
  std::vector<Schema::Column> projections_;
};

} // namespace MzPeak::Util
