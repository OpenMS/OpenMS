/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <algorithm>
#include <arrow/array.h>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "mzpeak/data/array_index.h"
#include "mzpeak/exception.h"
#include "mzpeak/schema/cv.h"
#include "mzpeak/util/slice.h"

// clang-format off
/**
 * Coordinate grid encoding (`MS:1003826`).
 *
 * A grid-encoded chunk stores no coordinate values.  Its `chunk_values` list is
 * null (or empty) and the coordinates are integer INDICES into a model whose
 * parameters ride in the same row, in a struct column
 * `<array>_grid { grid_type, parameters, indices }` that the array index
 * registers with `buffer_format = chunk_transform`, `transform = MS:1003826`.
 * The reader evaluates the model at each index.
 *
 * Grid encoding is chosen PER CHUNK.  The specification tells a writer that a
 * model whose error exceeds its threshold "SHOULD fall back to use a different
 * encoding", so one dimension may hold grid rows beside delta, Numpress or
 * uncompressed ones.  On the main axis each row's `chunk_encoding` says which;
 * on a secondary axis a row is grid when its grid struct is non-null.
 *
 * MAIN-axis indices are delta-coded -- `[first, deltas...]`, the start point
 * included; a secondary axis stores them as they are.
 *
 * Models and their parameter order:
 *
 *   MS:1003824  linear       [intercept, slope, scale=1]    (i*slope + intercept) / scale
 *   MS:1003825  square root  [intercept, slope, scale=1]    (i*slope + intercept)^2 / scale
 *   MS:9999002  timsTOF m/z  [C0, beta, C2, C3, C4, timebase, delay]
 *                            t = fma(i, timebase, delay);
 *                            solve t = C0 + beta*u + C2*u^2 (+ C3*u^3) for u;
 *                            m/z = u^2 - C4
 *   MS:9999001  TIMS 1/K0    [C6, C7, offset, slope]        1 / (C6 + C7/(offset + slope*i))
 *
 * The `MS:9999xxx` accessions are the reference's PLACEHOLDERS for the Bruker
 * models until PSI assigns terms; they are matched literally.
 *
 * WHICH ARITHMETIC.  There is no single reference to agree with bit for bit:
 * the reference grid codec (mzpeak_prototyping `grid.rs`) computes the flight
 * time with a fused multiply-add, while mzdata 0.66.7's scalar
 * `MzCalibrationModel2::convert_f64` multiplies and adds separately, and the
 * two differ in the last bit on roughly a third of timsTOF coordinates.  This
 * header follows `grid.rs`, because that is what WRITES grid-encoded archives:
 * a chunk's recorded bounds are its evaluation, so matching it is what makes
 * those bounds agree.  Every fusion is explicit (std::fma) and implicit
 * contraction is switched off per function below, so the result does not
 * depend on the compiler's default -- GCC contracts by default even in ISO
 * mode, and would otherwise fuse multiply-adds the reference keeps separate.
 */
// clang-format on

#if defined(__clang__)
#define MZPEAK_GRID_NO_CONTRACT _Pragma("clang fp contract(off)")
#else
#define MZPEAK_GRID_NO_CONTRACT
#endif

namespace MzPeak::Data::Transformer::Grid {

/// The chunk encoding this whole header exists for.
inline const Schema::CV& encoding_cv()
{
  static const Schema::CV cv("MS", "1003826");
  return cv;
}

/// Is @p entry the `<array>_grid` column of a grid-encoded dimension?
inline bool is_grid_encoded(const ArrayIndex::Entry& entry)
{
  return entry.buffer_format == Schema::BufferFormat::ChunkTransform &&
         entry.transform.has_value() && entry.transform->to_cv() == encoding_cv();
}

// GCC honours neither `STDC FP_CONTRACT` nor Clang's pragma, but it does honour
// an optimisation pragma scoped with push/pop, and that is all this needs: the
// model functions below, and nothing parsed after them.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("fp-contract=off")
#endif

/// Bruker timsTOF m/z at digitizer index @p i: the reference `grid.rs`
/// `timstof_mz`, operation for operation.
inline double timstof_mz(const std::vector<double>& p, double i)
{
  MZPEAK_GRID_NO_CONTRACT
  const double c0 = p[0], beta = p[1], c2 = p[2], c3 = p[3], c4 = p[4];
  const double timebase = p[5], delay = p[6];

  const double tof = std::fma(i, timebase, delay);
  const double s0 = (tof - c0) / beta;

  double refined = s0;
  if (c3 != 0.0) {
    // Cubic: Newton from the linear estimate, at most eight steps.
    double s = s0;
    const double c2_2 = c2 * 2.0, c3_3 = c3 * 3.0;
    for (int step = 0; step < 8; ++step) {
      const double s2 = s * s;
      const double f = std::fma(s, beta, c0) + s2 * c2 + std::fma(s2 * s, c3, -tof);
      const double deriv = std::fma(s, c2_2, beta) + s2 * c3_3;
      if (deriv == 0.0) break;
      const double delta = f / deriv;
      s -= delta;
      if (std::fabs(delta) < 1e-12) break;
    }
    refined = s;
  } else if (c2 != 0.0) {
    // Quadratic, closed form, the numerically stable root.
    const double d = beta * beta - 4.0 * c2 * (c0 - tof);
    refined = (d < 0.0) ? s0 : (c0 - tof) / (-0.5 * (beta + std::sqrt(d)));
  }
  return refined * refined - c4;
}

/// An order-of-magnitude fence around physically possible inverse reduced ion
/// mobility on a TIMS analyser, in Vs/cm^2.  Real acquisitions sit around
/// 0.6-1.9.  Not a tolerance: it exists to catch a mis-read parameter pair,
/// which lands near 45.
inline constexpr double kMobilityFloor = 0.1;
inline constexpr double kMobilityCeiling = 10.0;

/// Bruker TIMS 1/K0 at scan @p i: mzdata's `TimsCalibrationModel2::convert`.
///
/// `MS:9999001` carries `[c6, c7, offset, slope]` -- the order mzdata's
/// `as_param` writes and the reference's `from_param` reads.  Until upstream's
/// `eb08ba0` the reference WRITER emitted `[c6, c7, slope, offset]`, so writer
/// and reader disagreed, and nothing in such an archive says which it holds.
///
/// The pair is ordered by magnitude.  With `slope = (C3 - C2)/C1` and
/// `offset = C2 - slope*(C4 + C0)`, a real Bruker TIMS ramp -- decreasing
/// voltage (C2 > C3 >= 0), at least one scan (C1 >= 1), non-negative
/// C4 + C0 -- has |slope| <= C2 <= offset, so the larger magnitude is the
/// offset.  Those are PHYSICAL assumptions, not algebra: an increasing ramp can
/// violate them, and no acquisition seen uses one.  Because this is an
/// inference, every value it produces is checked against the fence above
/// rather than trusted; a wrong reading refuses the file instead of returning
/// plausible-looking numbers.  Equal magnitudes keep the declared order.
inline double timstof_mobility(const std::vector<double>& p, double i)
{
  MZPEAK_GRID_NO_CONTRACT
  const double c6 = p[0], c7 = p[1];
  const bool declared = std::fabs(p[2]) >= std::fabs(p[3]);
  const double offset = declared ? p[2] : p[3];
  const double slope = declared ? p[3] : p[2];
  return 1.0 / (c6 + c7 / (offset + slope * i));
}

/// The value of grid index @p idx under @p grid_type, or nullopt for a model
/// this reader does not know or a parameter list of the wrong length.
inline std::optional<double>
value_at(std::string_view grid_type, const std::vector<double>& p, std::uint32_t idx)
{
  MZPEAK_GRID_NO_CONTRACT
  const double i = static_cast<double>(idx);
  const std::size_t n = p.size();

  if (grid_type == "MS:1003824" && n >= 2 && n <= 3) {
    return (i * p[1] + p[0]) / (n == 3 ? p[2] : 1.0);
  }
  if (grid_type == "MS:1003825" && n >= 2 && n <= 3) {
    const double u = i * p[1] + p[0];
    return (u * u) / (n == 3 ? p[2] : 1.0);
  }
  if (grid_type == "MS:9999002" && n == 7) return timstof_mz(p, i);
  if (grid_type == "MS:9999001" && n == 4) return timstof_mobility(p, i);
  return std::nullopt;
}

#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif

namespace detail {

/// One string, tolerating both Arrow widths.
inline std::optional<std::string> string_at(const arrow::Array& a, int64_t i)
{
  if (a.IsNull(i)) return std::nullopt;
  if (a.type_id() == arrow::Type::LARGE_STRING) {
    return static_cast<const arrow::LargeStringArray&>(a).GetString(i);
  }
  if (a.type_id() == arrow::Type::STRING) {
    return static_cast<const arrow::StringArray&>(a).GetString(i);
  }
  return std::nullopt;
}

/// The values of row @p i of a list column, or null for a null row.
///
/// Accepts exactly the five list layouts `Util::Decoders::visit` does, so a
/// plain row decodes here whenever it would decode in a dimension with no grid
/// rows.  Anything else is refused: returning null for an unrecognised layout
/// used to hand a null array to the decoder, which dereferenced it.
inline std::shared_ptr<arrow::Array> list_row(const arrow::Array& a, int64_t i)
{
  if (!a.IsValid(i)) return nullptr;
  switch (a.type_id()) {
  case arrow::Type::LIST:
    return static_cast<const arrow::ListArray&>(a).value_slice(i);
  case arrow::Type::FIXED_SIZE_LIST:
    return static_cast<const arrow::FixedSizeListArray&>(a).value_slice(i);
  case arrow::Type::LARGE_LIST:
    return static_cast<const arrow::LargeListArray&>(a).value_slice(i);
  case arrow::Type::LIST_VIEW:
    return static_cast<const arrow::ListViewArray&>(a).value_slice(i);
  case arrow::Type::LARGE_LIST_VIEW:
    return static_cast<const arrow::LargeListViewArray&>(a).value_slice(i);
  default:
    throw TypeError("expected an arrow list array but found: " + a.type()->name());
  }
}

/// Read a numeric Arrow array as doubles (parameters are float64 by
/// specification; float32 is accepted rather than refused).
inline bool doubles_of(const arrow::Array& a, std::vector<double>& out)
{
  out.clear();
  out.reserve(static_cast<std::size_t>(a.length()));
  switch (a.type_id()) {
  case arrow::Type::DOUBLE:
    for (int64_t k = 0; k < a.length(); ++k) {
      out.push_back(static_cast<const arrow::DoubleArray&>(a).Value(k));
    }
    return true;
  case arrow::Type::FLOAT:
    for (int64_t k = 0; k < a.length(); ++k) {
      out.push_back(static_cast<const arrow::FloatArray&>(a).Value(k));
    }
    return true;
  default:
    return false;
  }
}

/// Read an integer Arrow array as grid indices.  The specification says
/// uint32; the narrower unsigned widths and int32 are accepted as the
/// reference accepts them.
inline bool indices_of(const arrow::Array& a, std::vector<std::uint32_t>& out)
{
  out.clear();
  out.reserve(static_cast<std::size_t>(a.length()));
  auto take = [&](const auto& typed) {
    for (int64_t k = 0; k < typed.length(); ++k) {
      out.push_back(static_cast<std::uint32_t>(typed.Value(k)));
    }
  };
  switch (a.type_id()) {
  case arrow::Type::UINT32:
    take(static_cast<const arrow::UInt32Array&>(a));
    return true;
  case arrow::Type::UINT16:
    take(static_cast<const arrow::UInt16Array&>(a));
    return true;
  case arrow::Type::UINT8:
    take(static_cast<const arrow::UInt8Array&>(a));
    return true;
  case arrow::Type::INT32:
    take(static_cast<const arrow::Int32Array&>(a));
    return true;
  default:
    return false;
  }
}

/// Flatten a per-row bound column across record batches; null stays null.
inline std::vector<std::optional<double>> flatten(const Util::Slice::Raw& raw)
{
  std::vector<std::optional<double>> out;
  for (const auto& chunk : raw) {
    const auto& a = *chunk;
    for (int64_t i = 0; i < a.length(); ++i) {
      if (a.IsNull(i)) {
        out.emplace_back();
      } else if (a.type_id() == arrow::Type::DOUBLE) {
        out.emplace_back(static_cast<const arrow::DoubleArray&>(a).Value(i));
      } else if (a.type_id() == arrow::Type::FLOAT) {
        out.emplace_back(static_cast<const arrow::FloatArray&>(a).Value(i));
      } else {
        throw InvalidFormatError("grid chunk bound column is not floating point");
      }
    }
  }
  return out;
}

/// Element @p i of @p v, or nullopt when the column was absent or shorter.
/// The bound columns are independent and either may be missing; indexing one
/// by the other's length is how a file without `chunk_end` read out of bounds.
inline std::optional<double> at(const std::vector<std::optional<double>>& v,
                                std::size_t i)
{
  return i < v.size() ? v[i] : std::nullopt;
}

} // namespace detail

/// One grid row, as read from a `<array>_grid` struct.
struct Row {
  std::string type;
  std::vector<double> parameters;
  std::vector<std::uint32_t> indices;
};

/**
 * Read row @p r of a grid struct batch into @p out.
 *
 * @returns false when the row is null, i.e. the chunk has no grid data.
 * @throws InvalidFormatError when the row is present but malformed.
 */
inline bool
read_row(const arrow::StructArray& rows, int64_t r, std::string_view dim, Row& out)
{
  if (rows.IsNull(r)) return false;

  auto grid_type = rows.GetFieldByName("grid_type");
  auto params = rows.GetFieldByName("parameters");
  auto idx = rows.GetFieldByName("indices");
  if (!grid_type || !params || !idx) {
    throw InvalidFormatError("grid column of " + std::string(dim) +
                             " lacks grid_type, parameters or indices");
  }

  const std::optional<std::string> type = detail::string_at(*grid_type, r);
  const std::shared_ptr<arrow::Array> p_row = detail::list_row(*params, r);
  const std::shared_ptr<arrow::Array> i_row = detail::list_row(*idx, r);
  if (!type || !p_row || !i_row || !detail::doubles_of(*p_row, out.parameters) ||
      !detail::indices_of(*i_row, out.indices)) {
    throw InvalidFormatError("malformed grid row in " + std::string(dim));
  }
  out.type = *type;
  return true;
}

/**
 * Evaluate a grid row, appending one coordinate per index to @p out.
 *
 * On return `row.indices` holds ABSOLUTE indices (the delta-coded main axis is
 * accumulated in place), which is what `step_at` needs.
 *
 * @throws InvalidFormatError for a model this reader does not know, for any
 *         non-finite coordinate, and for a `MS:9999001` mobility outside the
 *         physical fence.  A NaN would otherwise sail through every later
 *         check, since every comparison with it is false.
 */
inline void
evaluate(Row& row, bool delta_coded, std::string_view dim, std::vector<double>& out)
{
  if (delta_coded) {
    std::uint32_t acc = 0;
    for (auto& v : row.indices) {
      acc += v; // unsigned wrap-around, as the reference accumulates
      v = acc;
    }
  }

  const bool mobility = (row.type == "MS:9999001");
  for (const std::uint32_t k : row.indices) {
    const std::optional<double> value = value_at(row.type, row.parameters, k);
    if (!value) {
      throw InvalidFormatError("unknown grid model " + row.type + " with " +
                               std::to_string(row.parameters.size()) +
                               " parameters in " + std::string(dim));
    }
    if (!std::isfinite(*value)) {
      throw InvalidFormatError("grid model " + row.type + " gives a non-finite " +
                               "value at index " + std::to_string(k) + " in " +
                               std::string(dim));
    }
    if (mobility && !(*value > kMobilityFloor && *value < kMobilityCeiling)) {
      throw InvalidFormatError(
          "grid decode of " + std::string(dim) + " gives " + std::to_string(*value) +
          " Vs/cm^2 at scan " + std::to_string(k) +
          ", outside the physically possible range for inverse reduced ion "
          "mobility; the model's parameters cannot be read as written");
    }
    out.push_back(*value);
  }
}

/// How far a coordinate snapped to grid index @p idx can lie from the model's
/// value there: half the larger of the two adjacent grid steps.
///
/// A coordinate snaps to its NEAREST grid point, so it is at most half a step
/// away -- and a transcription error that shifts the model by one whole index,
/// which is exactly one step, is caught.  At one step it was not: a flight time
/// off by one digitizer tick passed every endpoint of the fixture.  Both
/// neighbours are consulted because the steps of a non-linear grid differ on
/// either side, and only the ones that exist: `idx + 1` on the largest index
/// used to wrap to zero and produce a "step" spanning the whole grid.
inline double snap_tolerance(const Row& row, std::uint32_t idx)
{
  const auto here = value_at(row.type, row.parameters, idx);
  if (!here) return 0.0;
  double step = 0.0;
  const auto widen = [&](std::uint32_t other) {
    if (const auto there = value_at(row.type, row.parameters, other)) {
      const double d = std::fabs(*there - *here);
      if (std::isfinite(d)) step = std::max(step, d);
    }
  };
  if (idx < std::numeric_limits<std::uint32_t>::max()) widen(idx + 1);
  if (idx > 0) widen(idx - 1);
  return 0.5 * step;
}

/**
 * Check a grid row's first and last coordinate against the bounds its chunk
 * recorded.
 *
 * The bounds are an oracle for free -- a writer that records its own
 * evaluation hands the reader a second opinion on every row, and a model
 * transcribed with one wrong operation still produces entirely plausible
 * numbers.  But the specification does not say bounds ARE evaluations: grid
 * encoding is "likely to be a lossy transformation" and a writer may record
 * the original coordinate instead.  Snapping to the nearest grid point moves a
 * coordinate by at most half a step, so that is the tolerance (see
 * snap_tolerance).  It admits every correctly snapped lossy grid and refuses a
 * model shifted by a whole index or more.  It does NOT catch an error smaller
 * than half a step, which only an exact reference can: exactness is enforced
 * by the tests in test/grid_encoding_test.cpp, not at read time.
 *
 * A recorded `chunk_end` of exactly 0.0 below a positive start, on a chunk of
 * ONE point, is treated as unrecorded.  The reference writer read a chunk's
 * start and end from a single iterator, so a one-value chunk had its value
 * consumed for the start and its end fell through to 0.0; its own
 * `chunk_series.rs` documents the fix.  In the conversion this was found on,
 * the two single-point chunks of 295 are exactly the two with a zero end.
 * The condition is kept to that shape: a multi-point chunk with a zero end is
 * a contradiction and is refused.
 */
inline void check_bounds(double first,
                         double last,
                         std::size_t points,
                         std::optional<double> start,
                         std::optional<double> end,
                         double tolerance_first,
                         double tolerance_last,
                         std::string_view dim,
                         std::size_t chunk)
{
  constexpr double kRelativeFloor = 1e-9;

  if (points == 1 && end && *end == 0.0 && start && *start > 0.0) end.reset();

  const auto check = [&](const std::optional<double>& expected, double got,
                         double snap, const char* which) {
    if (!expected) return;
    // A recorded bound of inf would make the tolerance below infinite and
    // accept any value at all, so a non-finite bound is refused outright.
    if (!std::isfinite(*expected)) {
      throw InvalidFormatError(std::string(which) + " of chunk " +
                               std::to_string(chunk) + " in " + std::string(dim) +
                               " is not finite");
    }
    const double tolerance =
        std::max(kRelativeFloor * std::max(std::fabs(*expected), 1.0), snap);
    // Written so a NaN on either side FAILS the check rather than passing it.
    if (!(std::fabs(got - *expected) <= tolerance)) {
      throw InvalidFormatError(
          "grid decode of " + std::string(dim) + " disagrees with " + which +
          " on chunk " + std::to_string(chunk) + ": model gives " +
          std::to_string(got) + ", file records " + std::to_string(*expected));
    }
  };
  check(start, first, tolerance_first, "chunk_start");
  check(end, last, tolerance_last, "chunk_end");
}

} // namespace MzPeak::Data::Transformer::Grid

#undef MZPEAK_GRID_NO_CONTRACT
