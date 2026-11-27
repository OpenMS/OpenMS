/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "mzpeak/spectrum_metadata.h"

namespace MzPeak {

/**
 * Per-chromatogram descriptive metadata.
 *
 * Read from the `chromatogram` struct column of a nested metadata table, or
 * from the top-level columns of a split-layout primary file.  The precursor and
 * selected-ion facets reuse the spectrum types: on disk they are the same
 * structures, joined the same way.
 *
 * @note There is deliberately no `products` member.  An SRM/MRM transition's Q3
 * isolation window belongs in a product facet, and the specification defines
 * one — but no reference writer emits it and the reference reader reaches an
 * unimplemented branch when it meets one.  Exposing an always-empty `products`
 * would suggest a chromatogram simply had no product, which for an SRM trace is
 * false.  Callers that need Q3 must be told it is unavailable, not handed an
 * empty vector; see @ref has_unreadable_product.
 */
struct ChromatogramMetadata final {
  /// `chromatogram.index` (uint64).
  uint64_t index = 0;

  /// `chromatogram.id` — the native chromatogram identifier
  /// (e.g. "TIC", "SRM SIC 500.2,184.1").
  std::string id;

  /// `MS_1000626_chromatogram_type` (CURIE, e.g. "MS:1000235" total ion
  /// current chromatogram, "MS:1001473" selected reaction monitoring
  /// chromatogram).  Empty string when null in the source.
  ///
  /// @note mzdata's `ChromatogramType::to_curie()` emits MS:1000472/MS:1000473
  /// for SIM/SRM while its own reader expects MS:1001472/MS:1001473.  Values
  /// are passed through unaltered; a caller comparing against one convention
  /// should accept both.
  std::string chromatogram_type;

  /// `MS_1000465_scan_polarity` (+1 positive, -1 negative).
  ///
  /// @note The reference writer emits 0 for "unknown" where the specification
  /// calls for null, so a 0 here means the source could not say.
  std::optional<int> polarity;

  /// `MS_1003060_number_of_data_points`.
  std::optional<uint64_t> number_of_data_points;

  /// `data_processing_ref` — reference to the data-processing chain applied to
  /// this chromatogram.  Empty string when null.
  std::string data_processing_ref;

  /// Per-chromatogram CV parameters (`chromatogram.parameters`).
  std::vector<CvParam> parameters;

  /// Precursor list.  For an SRM/MRM trace this carries Q1 (the selected ion
  /// and isolation window); Q3 lives in the product facet and is NOT available
  /// — see the note on this struct.
  std::vector<PrecursorInfo> precursors;

  /**
   * `true` when this chromatogram's type implies a product (Q3) selection that
   * the format cannot currently deliver — the selected/multiple REACTION
   * monitoring types.
   *
   * Selected ION monitoring (MS:1001472) is deliberately NOT flagged: it names
   * a precursor and no product, so reporting a lost product for it would be an
   * invention.  Its precursor is reachable through @ref precursors when the
   * file carries the facet.
   *
   * A caller reconstructing a transition MUST check this: the alternative is
   * silently treating an SRM trace as though it had no product, which produces
   * a plausible and wrong transition.
   */
  bool has_unreadable_product = false;
};

} // namespace MzPeak
