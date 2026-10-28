/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <arrow/array.h>
#include <arrow/builder.h>
#include <map>
#include <memory>
#include <ranges>
#include <string>
#include <vector>

#include "mzpeak/data/array_index.h"
#include "mzpeak/data/null_marking.h"
#include "mzpeak/data/signals.h"
#include "mzpeak/data/transformer/grid.h"
#include "mzpeak/data/transformer/primary.h"
#include "mzpeak/data/transformer/secondary.h"
#include "mzpeak/exception.h"
#include "mzpeak/schema/psi/data_type.h"
#include "mzpeak/util/decoders.h"
#include "mzpeak/util/slice.h"
#include "mzpeak/util/types.h"

namespace MzPeak::Data::Encoding {

/**
 * Decode mzPeak signal data encoding (point and chunk).
 *
 * The template type T should match the type for the main axis.
 */
template <typename T> class Decoder {
public:
  /// Constructor.
  Decoder(std::shared_ptr<Signals> signals,
          std::shared_ptr<Util::Slice> slice,
          const Util::DeltaEstimator<T>& estimator)
      : signals_(std::move(signals))
      , slice_(std::move(slice))
      , delta_estimator_(estimator)
  {
  }

  /**
   * Decode a float or double.
   */
  template <typename V>
  void decimal(const ArrayIndex::Dimension&, std::vector<V>&) const;

  /**
   * Decode a 32- or 64-bit integer.
   */
  template <typename V>
  void integer(const ArrayIndex::Dimension&, std::vector<V>&) const;

  /// The unit the values of @p dim were decoded in, or empty when the dimension
  /// declared none or drew rows from columns in two different units.  Load
  /// bearing: a time array in minutes vs seconds is a 60x error, and has_uv
  /// stores intensities in counts vs absorbance.  Only valid after decode.
  const std::string& unit_of(const ArrayIndex::Dimension& dim) const
  {
    static const std::string none;
    auto it = units_.find(dim.name);
    return it == units_.end() ? none : it->second;
  }

private:
  template <typename V>
  void decode(const ArrayIndex::Dimension&, std::vector<V>&) const;

  template <typename N, typename V>
  void decode_with_nulls(const ArrayIndex::Dimension&,
                         const N& null_decoder,
                         std::vector<V>&) const;

  template <Util::Type From, typename V>
  void remap(const ArrayIndex::Dimension& dim, std::vector<V>& v) const;

  /// Coordinate grid encoding (MS:1003826).  A grid row stores a STRUCT of
  /// model parameters and integer indices, not a list of values, so a dimension
  /// holding any grid rows cannot go through the list-shaped pipeline below and
  /// is decoded here -- all of its rows, grid and plain alike, because a
  /// consuming read leaves the columns to exactly one reader.
  template <typename N, typename V>
  void grid(const ArrayIndex::Dimension& dim,
            const N& null_decoder,
            std::vector<V>& v) const;

  void note_unit(const ArrayIndex::Dimension& dim, const std::string& unit) const
  {
    auto [it, inserted] = units_.try_emplace(dim.name, unit);
    if (!inserted && it->second != unit) it->second.clear();
  }

  template <typename V>
  using builder_for = arrow::NumericBuilder<
      typename Util::type_traits<Util::enum_type_v<V>>::array_type::TypeClass>;

  template <typename V>
  void finish(builder_for<V>& builder,
              const ArrayIndex::Dimension& dim,
              std::vector<V>& out) const;

  template <typename V>
  void reconstruct(const std::shared_ptr<arrow::Array>& src,
                   bool needs_delta_model,
                   std::vector<V>& out) const;

  /// One logical array spread across several complementary point columns
  /// (has_uv).  Upstream's single values_entry() path cannot express this.
  template <typename V>
  void coalesced_point(const ArrayIndex::Dimension&, std::vector<V>&) const;

  std::shared_ptr<Signals> signals_;
  std::shared_ptr<Util::Slice> slice_;
  Util::DeltaEstimator<T> delta_estimator_;
  mutable std::map<std::string, std::string> units_;
};

/******************************************************************************/
template <typename T>
template <typename V>
void Decoder<T>::decimal(const ArrayIndex::Dimension& dim, std::vector<V>& v) const
{
  Util::Type type = dim.type_or_throw();

  if (type == Util::Type::Float32) {
    remap<Util::Type::Float32>(dim, v);
  } else if (type == Util::Type::Float64) {
    remap<Util::Type::Float64>(dim, v);
  } else {
    Util::lift_type(type, []<Util::Type X> {
      std::string msg("Expected float or double but got: ");
      msg += Util::type_traits<X>::name;
      throw(TypeError(msg));
    });
  }
}

/******************************************************************************/
template <typename T>
template <typename V>
void Decoder<T>::integer(const ArrayIndex::Dimension& dim, std::vector<V>& v) const
{
  Util::Type type = dim.type_or_throw();

  if (type == Util::Type::Int32) {
    remap<Util::Type::Int32>(dim, v);
  } else if (type == Util::Type::Int64) {
    remap<Util::Type::Int64>(dim, v);
  } else {
    Util::lift_type(type, []<Util::Type X> {
      std::string msg("Expected int32 or int64 but got: ");
      msg += Util::type_traits<X>::name;
      throw(TypeError(msg));
    });
  }
}

/******************************************************************************/
template <typename T>
template <Util::Type From, typename V>
void Decoder<T>::remap(const ArrayIndex::Dimension& dim, std::vector<V>& v) const
{
  using F = Util::type_traits<From>::value_type;

  if constexpr (std::is_same_v<F, V>) {
    decode<V>(dim, v);
  } else {
    std::vector<F> tmp;
    decode<F>(dim, tmp);
    v.insert(v.end(), tmp.begin(), tmp.end());
  }
}

/******************************************************************************/
template <typename T>
template <typename V>
void Decoder<T>::decode(const ArrayIndex::Dimension& dim, std::vector<V>& v) const
{
  if (signals_->array_index()->layout() == ArrayIndex::Layout::Unknown) {
    throw UnknownLayoutError("cannot decode dimension, unknown layout: " + dim.name);
  }

  // Coalesced point: one logical array over several complementary point columns
  // (has_uv).  Upstream's dispatch decodes a single values_entry() and cannot
  // express this, so it is handled before the ordinary single-entry cases.
  const auto& entries = dim.entries;
  if (entries.size() > 1 && std::ranges::all_of(entries, [](const auto& e) {
        return e.buffer_format == Schema::BufferFormat::Point;
      })) {
    coalesced_point<V>(dim, v);
    return;
  }

  note_unit(dim, dim.values_entry().unit);

  switch (signals_->array_index()->layout()) {
  case ArrayIndex::Layout::Point:
  case ArrayIndex::Layout::Chunked:
    if (dim.needs_delta_model()) {
      using N = NullMarking::Decoder<V, T>;
      decode_with_nulls<N, V>(dim, N{delta_estimator_}, v);
    } else {
      using N = Util::Decoders::NullToZero<V>;
      decode_with_nulls<N, V>(dim, N{}, v);
    }
    break;

  case ArrayIndex::Layout::Unknown:
    throw UnknownLayoutError("cannot decode dimension: " + dim.name);
  }
}

/******************************************************************************/
template <typename T>
template <typename N, typename V>
void Decoder<T>::decode_with_nulls(const ArrayIndex::Dimension& dim,
                                   const N& null_decoder,
                                   std::vector<V>& v) const
{
  // Before `values_entry()`: on a secondary axis that mixes plain and grid rows
  // it names the plain column, and dispatching on it alone dropped the grid
  // rows without a word.
  if (std::ranges::any_of(dim.entries, &Transformer::Grid::is_grid_encoded)) {
    grid<N, V>(dim, null_decoder, v);
    return;
  }

  const auto& primary_entry = dim.values_entry();
  auto col = signals_->column(primary_entry);

  if (!col.has_value()) {
    throw InvalidFormatError("unable to decode dimension, not in schema: " +
                             dim.name);
  } else if (!slice_->has_column(col.value())) {
    return; // No data to decode so we can exit early.
  }

  auto go = [&](auto&& decoder) -> void { slice_->array(col.value(), v, decoder); };

  if (primary_entry.buffer_format == Schema::BufferFormat::Point) {
    auto decoder = Util::Decoders::Scalar<V, std::vector<V>, N>(null_decoder);
    go(decoder);
  } else {
    if (dim.is_main_axis()) {
      using Transformer = Transformer::Primary::Decoder<V>;
      Transformer transformer(signals_, slice_, dim);
      auto decoder = Util::Decoders::Flattened<V, std::vector<V>, N, Transformer>(
          null_decoder, std::move(transformer));
      go(decoder);
    } else {
      using Transformer = Transformer::Secondary::Decoder<V>;
      Transformer transformer(dim);
      auto decoder = Util::Decoders::Flattened<V, std::vector<V>, N, Transformer>(
          null_decoder, std::move(transformer));
      go(decoder);
    }
  }
}

/******************************************************************************/
template <typename T>
template <typename N, typename V>
void Decoder<T>::grid(const ArrayIndex::Dimension& dim,
                      const N& null_decoder,
                      std::vector<V>& v) const
{
  namespace Grid = Transformer::Grid;
  using BF = Schema::BufferFormat;
  using Encoding = Schema::PSI::ChunkEncoding;
  using Transform = Schema::PSI::Transform;

  const bool main_axis = dim.is_main_axis();
  const auto grid_entry = std::ranges::find_if(dim.entries, &Grid::is_grid_encoded);

  // The dimension's OTHER chunk_transform column, if it has one.  That is where
  // Numpress keeps its bytes: `small.numpress.mzpeak` carries them in
  // `mz_numpress_linear_bytes` with every `mz_chunk_values` null, so reading a
  // Numpress row from chunk_values finds nothing at all.
  const auto bytes_entry = std::ranges::find_if(dim.entries, [](const auto& e) {
    return e.buffer_format == BF::ChunkTransform && !Grid::is_grid_encoded(e);
  });
  const bool has_bytes = bytes_entry != dim.entries.end();

  if (!main_axis) {
    // A secondary Numpress transform also keeps its bytes in a column of its
    // own, and applying it row by row beside grid rows is not implemented;
    // refusing is better than returning a mobility array with holes in it.
    // The zero-intensity transforms are pass-throughs -- the null decoder does
    // their work -- so a plain column carrying one is decoded as usual.
    const auto plain = dim.entry_with(BF::ChunkSecondary);
    const auto t = plain && plain->transform ? plain->transform->type()
                                             : std::optional<Transform::Type>();
    const bool pass_through = !plain || !plain->transform ||
                              t == Transform::ZeroIntensityTrim ||
                              t == Transform::ZeroIntensityInterpolation;
    if (has_bytes || !pass_through) {
      throw InvalidFormatError(
          "dimension " + dim.name +
          " mixes grid rows with a Numpress-transformed "
          "secondary column, which this reader does not decode");
    }
  }

  // Every column is taken exactly once, here: a Slice read consumes it.
  const auto take = [&](const std::optional<Schema::Column>& c) {
    return (c && slice_->has_column(*c)) ? slice_->raw(*c)
                                         : std::shared_ptr<Util::Slice::Raw>();
  };
  const std::shared_ptr<Util::Slice::Raw> grid_raw =
      take(signals_->column(*grid_entry));
  const std::shared_ptr<Util::Slice::Raw> values_raw =
      take(signals_->column(dim, main_axis ? BF::ChunkValues : BF::ChunkSecondary));
  const std::shared_ptr<Util::Slice::Raw> bytes_raw =
      main_axis && has_bytes ? take(signals_->column(*bytes_entry))
                             : std::shared_ptr<Util::Slice::Raw>();
  if (!grid_raw && !values_raw && !bytes_raw) return;

  std::vector<std::string> encodings;
  std::vector<std::optional<double>> starts, ends;
  if (main_axis) {
    // Required, as the ordinary chunk decoder requires it.  Inferring a main
    // axis row's kind from its grid struct -- what a secondary axis has to do
    // -- would make a missing column an error only for files with plain rows.
    const auto enc = signals_->column(dim, BF::ChunkEncoding);
    if (!enc || !slice_->has_column(*enc)) {
      throw InvalidFormatError("while decoding " + dim.name +
                               " a needed column with buffer format "
                               "chunk_encoding was not found");
    }
    slice_->array(*enc, encodings, Util::Decoders::Scalar<std::string>());
    if (auto raw = take(signals_->column(dim, BF::ChunkStart))) {
      starts = Grid::detail::flatten(*raw);
    }
    if (auto raw = take(signals_->column(dim, BF::ChunkEnd))) {
      ends = Grid::detail::flatten(*raw);
    }
  }

  const auto encoding_of = [&](std::size_t flat) -> Encoding::Type {
    if (flat >= encodings.size()) {
      throw InvalidFormatError("while decoding " + dim.name + " chunk " +
                               std::to_string(flat) + " has no chunk_encoding");
    }
    const auto cv = Schema::CV::from_string(encodings[flat]);
    if (!cv) {
      throw InvalidFormatError("invalid chunk encoding CV: " + encodings[flat]);
    }
    const auto type = Encoding(*cv).type();
    if (!type) {
      throw InvalidFormatError("while decoding " + dim.name +
                               " unknown chunk encoding method: " + encodings[flat]);
    }
    return *type;
  };

  // The three columns come from the same record batches; say so if they don't.
  const std::shared_ptr<Util::Slice::Raw> sources[] = {grid_raw, values_raw,
                                                       bytes_raw};
  std::size_t batches = 0;
  for (const auto& src : sources) {
    if (src) batches = src->size();
  }
  for (const auto& src : sources) {
    if (src && src->size() != batches) {
      throw InvalidFormatError("the chunk columns of " + dim.name +
                               " are batched differently");
    }
  }

  Util::Decoders::Scalar<V, std::vector<V>, N> scalar(null_decoder);
  Grid::Row row;
  std::vector<double> coords;
  std::optional<double> previous_start, previous_end;
  std::size_t flat = 0;

  for (std::size_t b = 0; b < batches; ++b) {
    const arrow::StructArray* g = nullptr;
    if (grid_raw) {
      if ((*grid_raw)[b]->type_id() != arrow::Type::STRUCT) {
        throw InvalidFormatError("grid column of " + dim.name + " is not a struct");
      }
      g = static_cast<const arrow::StructArray*>((*grid_raw)[b].get());
    }
    const arrow::Array* values = values_raw ? (*values_raw)[b].get() : nullptr;
    const arrow::Array* bytes = bytes_raw ? (*bytes_raw)[b].get() : nullptr;

    const int64_t length = g        ? g->length()
                           : values ? values->length()
                                    : bytes->length();
    for (const arrow::Array* a :
         {static_cast<const arrow::Array*>(g), values, bytes}) {
      if (a && a->length() != length) {
        throw InvalidFormatError("the chunk columns of " + dim.name +
                                 " disagree on row count");
      }
    }

    for (int64_t r = 0; r < length; ++r, ++flat) {
      const std::size_t first = v.size();

      // On the main axis a row's chunk_encoding says what it is.  A secondary
      // axis has none of its own, so a row is grid when its grid struct is
      // present -- and a row carrying both has no basis for preferring one.
      std::optional<Encoding::Type> encoding;
      bool is_grid = false;
      if (main_axis) {
        encoding = encoding_of(flat);
        is_grid = *encoding == Encoding::Type::Grid;
      } else {
        const bool has_grid = g && g->IsValid(r);
        const auto plain = values ? Grid::detail::list_row(*values, r) : nullptr;
        if (has_grid && plain && plain->length() > 0) {
          throw InvalidFormatError("chunk " + std::to_string(flat) + " of " +
                                   dim.name +
                                   " carries both plain values and a grid");
        }
        is_grid = has_grid;
      }

      if (is_grid) {
        if (!g || !Grid::read_row(*g, r, dim.name, row)) {
          if (main_axis) {
            throw InvalidFormatError("chunk " + std::to_string(flat) + " of " +
                                     dim.name +
                                     " declares grid encoding but has no grid");
          }
          continue;
        }
        coords.clear();
        Grid::evaluate(row, main_axis, dim.name, coords);
        if (main_axis && !coords.empty()) {
          Grid::check_bounds(
              coords.front(), coords.back(), coords.size(),
              Grid::detail::at(starts, flat), Grid::detail::at(ends, flat),
              Grid::snap_tolerance(row, row.indices.front()),
              Grid::snap_tolerance(row, row.indices.back()), dim.name, flat);
        }
        for (const double c : coords)
          v.push_back(static_cast<V>(c));
      } else if (!main_axis) {
        if (values) {
          if (const auto src = Grid::detail::list_row(*values, r))
            scalar.decode(src, v);
        }
      } else {
        // Exactly what Transformer::Primary::Decoder::operator() and the
        // Flattened decoder do with a row -- from the column that encoding
        // keeps its data in -- so a plain row decodes here as it would in a
        // dimension with no grid rows at all.  A null row carries no data, as
        // the ordinary decoder's visit() treats it.
        const arrow::Array* source =
            *encoding == Encoding::Type::NumpressLinear ? bytes : values;
        const auto src = source ? Grid::detail::list_row(*source, r) : nullptr;
        if (src) {
          const auto start = Grid::detail::at(starts, flat);
          if (!start) {
            throw InvalidFormatError("chunk " + std::to_string(flat) + " of " +
                                     dim.name + " lacks its chunk_start");
          }
          switch (*encoding) {
          case Encoding::Type::NoCompression:
            v.push_back(static_cast<V>(*start));
            scalar.decode(src, v);
            break;
          case Encoding::Type::Delta:
            scalar.decode(Util::Algorithm::null_delta_decode<Util::enum_type_v<V>>(
                              static_cast<V>(*start), src),
                          v);
            break;
          case Encoding::Type::NumpressLinear: {
            const auto decoded = Util::Numpress::decode_linear_convert<V>(src);
            v.insert(v.end(), decoded->begin(), decoded->end());
            break;
          }
          case Encoding::Type::Grid:
            std::unreachable(); // routed above
          }
        }
      }

      // The structural checks Primary::Decoder makes on every chunk -- start
      // <= end, chunks ascending and non-overlapping -- applied to grid and
      // plain rows alike.  Primary skips a row with no end, which let a
      // descending chunk through after it; here the order of STARTS is always
      // enforced, and a missing or excused end is replaced by the row's last
      // decoded value.  A one-point chunk with a zero end is the writer bug
      // check_bounds() describes; start == end == 0 marks an empty chunk.
      if (main_axis) {
        const auto s = Grid::detail::at(starts, flat);
        auto e = Grid::detail::at(ends, flat);
        const std::size_t points = v.size() - first;
        if (points == 1 && e && *e == 0.0 && s && *s > 0.0) e.reset();
        const bool empty_marker = s && e && *s == 0.0 && *e == 0.0;
        if (s && !empty_marker) {
          if (e && !(*s <= *e)) {
            throw InvalidFormatError("chunk_start exceeds chunk_end (" + dim.name +
                                     ")");
          }
          if ((previous_start && *s < *previous_start) ||
              (previous_end && *s < *previous_end)) {
            throw InvalidFormatError("chunks overlap or are out of order (" +
                                     dim.name + ")");
          }
          previous_start = s;
          if (e) {
            previous_end = e;
          } else if (points > 0) {
            previous_end = static_cast<double>(v.back());
          }
        }
      }
    }
  }
}

/******************************************************************************/
template <typename T>
template <typename V>
void Decoder<T>::finish(builder_for<V>& builder,
                        const ArrayIndex::Dimension& dim,
                        std::vector<V>& out) const
{
  std::shared_ptr<arrow::Array> assembled;
  if (!builder.Finish(&assembled).ok()) {
    throw ParquetError("failed to assemble array for dimension: " + dim.name);
  }
  reconstruct<V>(assembled, dim.needs_delta_model(), out);
}

/******************************************************************************/
template <typename T>
template <typename V>
void Decoder<T>::reconstruct(const std::shared_ptr<arrow::Array>& src,
                             bool needs_delta_model,
                             std::vector<V>& out) const
{
  if (needs_delta_model) {
    using N = NullMarking::Decoder<V, T>;
    Util::Decoders::Scalar<V, std::vector<V>, N> decoder{N{delta_estimator_}};
    decoder.decode(src, out);
  } else {
    using N = Util::Decoders::NullToZero<V>;
    Util::Decoders::Scalar<V, std::vector<V>, N> decoder{N{}};
    decoder.decode(src, out);
  }
}

/******************************************************************************/
template <typename T>
template <typename V>
void Decoder<T>::coalesced_point(const ArrayIndex::Dimension& dim,
                                 std::vector<V>& v) const
{
  using array_type = Util::type_traits<Util::enum_type_v<V>>::array_type;

  std::vector<const ArrayIndex::Entry*> ordered;
  ordered.reserve(dim.entries.size());
  for (const auto& e : dim.entries) {
    if (e.buffer_priority) ordered.push_back(&e);
  }
  for (const auto& e : dim.entries) {
    if (!e.buffer_priority) ordered.push_back(&e);
  }

  std::vector<std::shared_ptr<Util::Slice::Raw>> columns;
  std::vector<std::string> column_units;
  for (const auto* e : ordered) {
    auto col = signals_->array_index()->entry_column(*signals_->groups(), *e);
    if (!col.has_value()) continue;
    if (auto raw = slice_->raw(col.value())) {
      columns.push_back(std::move(raw));
      column_units.push_back(e->unit);
    }
  }

  if (columns.empty()) {
    throw ParquetError("unable to decode dimension, not in schema: " + dim.name);
  }

  builder_for<V> builder;

  const std::size_t n_chunks = columns.front()->size();
  for (std::size_t c = 0; c < n_chunks; ++c) {
    auto first = std::static_pointer_cast<array_type>((*columns.front())[c]);
    for (int64_t r = 0; r < first->length(); ++r) {
      bool written = false;
      for (std::size_t ci = 0; ci < columns.size(); ++ci) {
        const auto& column = columns[ci];
        if (c >= column->size()) continue;
        auto arr = std::static_pointer_cast<array_type>((*column)[c]);
        if (r >= arr->length() || arr->IsNull(r)) continue;
        if (written) {
          throw ParquetError(
              "two columns of dimension '" + dim.name +
              "' carry a value on the same row; they may be in different units "
              "and there is no basis for preferring one");
        }
        (void)builder.Append(arr->Value(r));
        note_unit(dim, column_units[ci]);
        written = true;
      }
      if (!written) (void)builder.AppendNull();
    }
  }

  finish<V>(builder, dim, v);
}

} // namespace MzPeak::Data::Encoding
