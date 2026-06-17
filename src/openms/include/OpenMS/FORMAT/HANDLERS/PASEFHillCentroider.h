// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace OpenMS::Internal::PASEFHillCentroider
{
  /// Parameters for IM-axis hill-based centroiding of one TIMS-PASEF frame.
  /// Algorithm modeled on Biosaur2 (linkScanToHills_/splitHills_) but with
  /// IM scans within a frame playing the role of consecutive RT spectra.
  struct OPENMS_DLLAPI Params
  {
    /// m/z linking tolerance in ppm. A peak in the next IM scan extends a
    /// hill if its m/z is within this tolerance of the previous-scan peak.
    double mz_tol_ppm = 8.0;

    /// Hill-valley factor (hvf): both (left_max / valley) and
    /// (right_max / valley) must exceed this for a valley to split a hill.
    /// Smaller values split more aggressively. Biosaur2 default = 1.3.
    double valley_factor = 1.3;

    /// Minimum peaks per hill. Hills shorter than this are dropped.
    ///
    /// Default 1 is appropriate for TIMS-PASEF data: detector centroiding
    /// produces sparse per-scan peaks, and many real ions span only one IM
    /// scan in a frame. Empirically, on detector-centroided MS1 frames
    /// roughly 75% of peaks have no same-m/z partner in the previous IM
    /// scan within 100 ppm — the data is genuinely sparse along the IM
    /// axis, not noisy. Raising to 2 (or higher) is appropriate when the
    /// caller knows multi-scan hills are typical (e.g. RT-axis hill
    /// detection on dense LC-MS data, the original Biosaur2 use case).
    Size min_hill_length = 1;

    /// Discard input peaks below this intensity before linking.
    double min_intensity = 0.0;

    /// Tolerance for treating two flat-input IM values as belonging to the
    /// same IM scan. 0.0 means exact equality. TIMS calibration produces
    /// discrete IM values per scan index, so exact equality usually works;
    /// callers that have float-cast their IM values may want a small epsilon.
    double im_group_tolerance = 0.0;

    /// Maximum number of consecutive empty IM scans a hill is allowed to
    /// bridge while linking. 0 (default) = strict consecutive-scan linking
    /// (original Biosaur2 behavior). 1 = a single empty scan at the hill's
    /// m/z is tolerated; the hill's tail stays "alive" for one extra scan
    /// and can still be extended at scan+2. Useful for detector-centroided
    /// TIMS-PASEF where ions occasionally fail to register in one IM scan.
    /// Note: hill length counts only the scans where the ion was actually
    /// observed, not the bridged gap.
    ///
    /// REQUIRES @p scan_indices to be passed to centroidFrame(). Without it
    /// the helper cannot tell which IM scans are empty (vs. simply absent
    /// from the input), and silently treats this as 0 to avoid merging
    /// hills across unknown gaps.
    Size max_scan_gap = 0;
  };

  /// One centroided peak produced by hill detection.
  struct OPENMS_DLLAPI Centroid
  {
    double mz;          ///< intensity-weighted m/z over the hill
    double intensity;   ///< summed intensity across the hill
    double im_apex;     ///< IM value at the hill's intensity apex
    double im_lower;    ///< minimum IM value covered by the hill
    double im_upper;    ///< maximum IM value covered by the hill
    double mz_lower;    ///< minimum m/z value covered by the hill
    double mz_upper;    ///< maximum m/z value covered by the hill
    Size   length;      ///< number of IM scans contributing to the hill
  };

  /// Centroid one TIMS-PASEF frame by IM-axis hill detection.
  ///
  /// The input is a flat list of peaks from a single frame. The helper
  ///   1. groups peaks by IM scan (consecutive equal-im values),
  ///   2. links peaks across consecutive IM scans by m/z proximity (greedy
  ///      assignment in descending current-scan intensity order, mirroring
  ///      Biosaur2's linkScanToHills_),
  ///   3. splits hills at intensity valleys (port of Biosaur2's splitHills_),
  ///   4. drops hills shorter than @p p.min_hill_length, and
  ///   5. collapses each surviving hill to one Centroid.
  ///
  /// @param mz_values     Array of length @p n.
  /// @param intensities   Array of length @p n.
  /// @param im_values     Array of length @p n. Need not be sorted; the helper
  ///                      sorts internally. All peaks at one IM scan must
  ///                      share an im_value to within @p
  ///                      Params::im_group_tolerance (used only when
  ///                      @p scan_indices is null).
  /// @param n             Number of input peaks.
  /// @param p             Parameters.
  /// @param scan_indices  Optional. If provided, an array of length @p n
  ///                      giving the integer TIMS IM-scan index of each
  ///                      input peak. When provided, the linker iterates
  ///                      every scan index from min to max (including
  ///                      scans with no peaks) which makes @p
  ///                      Params::max_scan_gap meaningful — empty scans
  ///                      age active hill tails and a hill is closed once
  ///                      its tail's age exceeds max_scan_gap. When null,
  ///                      grouping falls back to im_values and gap
  ///                      detection is disabled.
  /// @return Centroids sorted by m/z ascending. Empty if @p n == 0 or all
  ///         peaks fail the filters.
  OPENMS_DLLAPI std::vector<Centroid>
  centroidFrame(const double* mz_values,
                const double* intensities,
                const double* im_values,
                std::size_t n,
                const Params& p,
                const std::uint32_t* scan_indices = nullptr);

} // namespace OpenMS::Internal::PASEFHillCentroider
