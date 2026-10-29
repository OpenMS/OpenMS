/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

namespace MzPeak {

/**
 * TOF -> m/z calibration for the Bruker TDF "ims-compact" layout.
 *
 * That layout stores NO m/z array.  A non-standard `tof` (Int32) column stands
 * in for it and m/z is reconstructed as `(a + b * tof)^2`, with the
 * coefficients declared in the archive index under `metadata.ims_calibration`.
 *
 * A reader that does not do this finds no m/z dimension at all, so its
 * m/z + intensity filter matches nothing and EVERY spectrum reads as empty —
 * a total, silent data loss that looks like an empty file rather than an
 * unsupported layout.
 *
 * @ref valid is false when the archive declares no calibration; in that case
 * `tof` is left alone rather than converted with meaningless coefficients.
 */
struct ImsCalibration final {
  /// Intercept of the calibration.
  double a = 0.0;

  /// Slope of the calibration.
  double b = 0.0;

  /// True once both coefficients have been read from the index.
  bool valid = false;

  /// Convert a TOF index to m/z.
  double mz(double tof) const
  {
    const double root = a + b * tof;
    return root * root;
  }
};

} // namespace MzPeak
