/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <arrow/array.h>

#include "mzpeak/data/array_index.h"
#include "mzpeak/schema/group.h"
#include "mzpeak/util/parquet.h"
#include "mzpeak/util/query.h"

namespace MzPeak::Data {

/**
 * Low-level access to a single data table in a Parquet file.
 */
class Signals {
public:
  /// Constructor.
  Signals(std::unique_ptr<Util::Parquet> parquet);

  /// Destructor.
  ~Signals();

  /**
   * Access the ArrayIndex for this data file.
   */
  const std::shared_ptr<ArrayIndex>& array_index() const;

  /**
   * Get the number of records in this data file.
   */
  std::size_t record_count() const;

  /**
   * Return a query builder for the index column.
   *
   * Useful for index-based queries such as `index().eq(n)`.
   */
  Util::Query::Builder index() const;

  /**
   * Execute a query, projecting the requested dimensions.
   */
  std::unique_ptr<Util::Slice> select(const std::vector<ArrayIndex::Dimension>&,
                                      const Util::Query&);

  /**
   * Low-level interface for accessing a column given its name.
   *
   * Useful if you need to manually construct queries.
   */
  std::optional<Schema::Column> column(std::string_view) const;

  /**
   * Low-level interface for accessing a column given an array index entry.
   */
  std::optional<Schema::Column> column(const ArrayIndex::Entry&) const;

  /**
   * Low-level interface for accessing a column given a dimension and
   * a buffer format.
   *
   * Returns the first matching column.
   */
  std::optional<Schema::Column> column(const ArrayIndex::Dimension&,
                                       Schema::BufferFormat) const;

  /**
   * Low-level interface for accessing the schema encoded as a map of
     Group objects.
   */
  const std::shared_ptr<Schema::GroupMap>& groups() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

} // namespace MzPeak::Data
