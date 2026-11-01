/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <arrow/array/array_primitive.h>
#include <arrow/io/api.h>
#include <arrow/result.h>
#include <arrow/table.h>
#include <arrow/util/config.h>
#include <memory>
#include <parquet/arrow/reader.h>
#include <vector>

#include "mzpeak/io/file.h"

namespace MzPeak::Util {

/**
 * ReadTable through whichever overload this Arrow has.
 *
 * The Result-returning forms arrived in Arrow 24 and deprecated the
 * out-parameter ones.  bioconda's OpenMS 3.5.0 pins Arrow 21, where only the
 * out-parameter forms exist -- and a consumer that links both this library and
 * that OpenMS must use ONE Arrow.  One shim here, used by the library and its
 * tests, rather than a version floor that rules the packaged OpenMS out.
 *
 * @param columns  leaf column indices to read, or null for every column.
 */
inline arrow::Result<std::shared_ptr<arrow::Table>>
read_table(parquet::arrow::FileReader& reader,
           const std::vector<int>* columns = nullptr)
{
#if ARROW_VERSION_MAJOR >= 24
  return columns ? reader.ReadTable(*columns) : reader.ReadTable();
#else
  std::shared_ptr<arrow::Table> table;
  const arrow::Status status =
      columns ? reader.ReadTable(*columns, &table) : reader.ReadTable(&table);
  if (!status.ok()) return status;
  return table;
#endif
}

/**
 * This is a low-level interface for accessing a Parquet file.
 */
class Arrow final {
public:
  using random_access_t = arrow::io::RandomAccessFile;

  /// Constructor.
  Arrow(std::unique_ptr<IO::File>);

  /// Destructor.
  ~Arrow();

  /**
   * Return an arrow I/O object that can be used to open a Parquet
   * file for reading.
   */
  std::shared_ptr<random_access_t> reader() const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

/// Cast to a typed numeric array (arrow::NumericArray<T>) by the array's type id. Those arrays are
/// templates that libarrow does not export, so a dynamic_cast to one fails across the library
/// boundary on macOS (libc++ compares type_info by address): the value is read as absent there.
template <typename ArrayT>
std::shared_ptr<ArrayT> numeric_array_cast(const std::shared_ptr<arrow::Array>& array)
{
  if (!array || array->type_id() != ArrayT::TypeClass::type_id) return nullptr;
  return std::static_pointer_cast<ArrayT>(array);
}

} // namespace MzPeak::Util
