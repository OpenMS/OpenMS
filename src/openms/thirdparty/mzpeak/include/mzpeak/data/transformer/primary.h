/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <utility>

#include "mzpeak/data/array_index.h"
#include "mzpeak/data/signals.h"
#include "mzpeak/schema/psi/chunk_encoding.h"
#include "mzpeak/util/algorithm.h"
#include "mzpeak/util/numpress.h"
#include "mzpeak/util/slice.h"
#include "mzpeak/util/types.h"

namespace MzPeak::Data::Transformer::Primary {

/**
 * Transform Arrow arrays that have been encoded with one of the chunk
 * encoding schemes on the primary axis.
 */
template <Util::supported_type T> class Decoder final {
public:
  /// The type of result this function object returns.
  using result_type = std::variant<std::shared_ptr<arrow::Array>,
                                   std::pair<T, std::shared_ptr<arrow::Array>>,
                                   std::shared_ptr<std::vector<T>>>;

  /// Constructor.
  Decoder(std::shared_ptr<Signals>,
          std::shared_ptr<Util::Slice>,
          const ArrayIndex::Dimension&);

  /**
   * Transform the given array.
   *
   * This array could come from the `chunk_values` column, or from the
   * `chunk_transform` column.
   */
  result_type operator()(int64_t, const std::shared_ptr<arrow::Array>&) const;

private:
  std::string dim_name_;
  std::vector<Schema::PSI::ChunkEncoding> chunk_encoding_;
  std::vector<T> chunk_start_;
};

/******************************************************************************/
template <Util::supported_type T>
Decoder<T>::Decoder(std::shared_ptr<Signals> signals,
                    std::shared_ptr<Util::Slice> slice,
                    const ArrayIndex::Dimension& dim)
    : dim_name_(dim.name)
    , chunk_encoding_()
    , chunk_start_()
{
  auto decode = [&]<typename U>(Schema::BufferFormat format,
                                std::vector<U>& dest) -> void {
    std::optional<Schema::Column> column = signals->column(dim, format);

    if (!column.has_value()) {
      std::string msg("while decoding " + dim.name);
      msg += " a needed column with buffer format ";
      msg += Schema::buffer_format_to_string(format);
      msg += " was not found";
      throw InvalidFormatError(msg);
    }

    slice->array(column.value(), dest, Util::Decoders::Scalar<U>());
  };

  // Decode the `chunk_encoding` column.
  std::vector<std::string> encodings;
  decode(Schema::BufferFormat::ChunkEncoding, encodings);
  chunk_encoding_.reserve(encodings.size());

  for (const auto& s : encodings) {
    std::optional<Schema::CV> cv = Schema::CV::from_string(s);

    if (!cv.has_value()) {
      throw InvalidFormatError("invalid chunk encoding CV: " + std::string(s));
    }

    chunk_encoding_.emplace_back(cv.value());
  }

  // Decode the `chunk_start` column.
  decode(Schema::BufferFormat::ChunkStart, chunk_start_);

  // Sanity checks:
  if (chunk_encoding_.size() != chunk_start_.size()) {
    std::string msg("while decoding " + dim.name);
    msg += " the chunk encoding and chunk start columns";
    msg += " have different lengths";
    throw InvalidFormatError(msg);
  }

  // Reject a file whose chunks are not ascending by chunk_start and
  // non-overlapping, or whose chunk_start exceeds its chunk_end.  The spec
  // requires this, and nothing else checks it, so a malformed file would
  // otherwise decode to plausible wrong values.  Done here, alongside the
  // chunk_start read, because Slice::raw() is a CONSUMING read -- validating in
  // a separate pass would erase chunk_start before this decode could read it.
  //
  // chunk_end is optional (some files omit it); its absence just skips the
  // check.  Read via the raw arrays so a null chunk_end -- the empty-chunk
  // marker -- is distinguishable from a genuine 0 rather than misread.
  std::optional<Schema::Column> end_col =
      signals->column(dim, Schema::BufferFormat::ChunkEnd);
  if (end_col.has_value() && slice->has_column(end_col.value())) {
    std::shared_ptr<Util::Slice::Raw> ends_raw = slice->raw(end_col.value());
    std::optional<T> previous_end;
    std::size_t flat = 0;
    if (ends_raw != nullptr) {
      for (const auto& chunk : *ends_raw) {
        auto ends = std::static_pointer_cast<
            typename Util::type_traits<Util::enum_type_v<T>>::array_type>(chunk);
        for (int64_t r = 0; r < ends->length(); ++r, ++flat) {
          if (ends->IsNull(r) || flat >= chunk_start_.size()) continue;
          const T start = chunk_start_[flat];
          const T end = ends->Value(r);
          // start == end == 0 is the empty-chunk sentinel.
          if (start == T{} && end == T{}) continue;
          if (!(start <= end)) {
            throw InvalidFormatError("chunk_start exceeds chunk_end (" + dim.name +
                                     ")");
          }
          if (previous_end.has_value() && start < *previous_end) {
            throw InvalidFormatError("chunks overlap or are out of order (" +
                                     dim.name + ")");
          }
          previous_end = end;
        }
      }
    }
  }
}

/******************************************************************************/
template <Util::supported_type T>
Decoder<T>::result_type
Decoder<T>::operator()(int64_t index, const std::shared_ptr<arrow::Array>& src) const
{
  if (index < 0 || static_cast<std::size_t>(index) >= chunk_encoding_.size()) {
    std::string msg("while decoding " + dim_name_);
    msg += " the chunk_values/chunk_transform column is out of bounds ";
    msg += std::to_string(index) + " >= " + std::to_string(chunk_encoding_.size());
    throw InvalidFormatError(msg);
  }

  std::optional<Schema::PSI::ChunkEncoding::Type> type =
      chunk_encoding_[index].type();

  if (!type.has_value()) {
    std::string msg("while decoding " + dim_name_);
    msg += " unknown chunk encoding method: ";
    msg += chunk_encoding_[index].to_cv().to_string();
    throw InvalidFormatError(msg);
  }

  switch (type.value()) {
  case Schema::PSI::ChunkEncoding::Type::NoCompression:
    return std::make_pair(chunk_start_[index], src);
  case Schema::PSI::ChunkEncoding::Type::Delta:
    return Util::Algorithm::null_delta_decode<Util::enum_type_v<T>>(
        chunk_start_[index], src);
  case Schema::PSI::ChunkEncoding::Type::NumpressLinear:
    return Util::Numpress::decode_linear_convert<T>(src);
  case Schema::PSI::ChunkEncoding::Type::Grid:
    // A grid-encoded dimension is routed to Transformer::Grid before this
    // decoder is ever built; reaching here means a chunk_values list was
    // populated on a dimension that declared no grid column to go with it.
    throw InvalidFormatError("while decoding " + dim_name_ +
                             " a chunk declares grid encoding but the "
                             "dimension has no grid column");
  }

  std::unreachable();
}

} // namespace MzPeak::Data::Transformer::Primary
