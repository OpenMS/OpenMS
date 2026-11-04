/*

This file is part of the mzpeak.h project.  It is subject to the
license specified in the LICENSE file which can be found in the
top-level directory of this repository.

*/

#pragma once

#include <memory>
#include <optional>
#include <vector>

namespace MzPeak::Util {
class Parquet;
}

namespace MzPeak::Metadata {

/**
 * Metadata from the spectrum table.
 */
class Spectrum final {
public:
  /// Constructor.
  Spectrum(std::unique_ptr<Util::Parquet>, uint64_t);

  /**
   * Return the spectrum level (MS:1000511).
   *
   * NOTE: If no level was found in the metadata this function will
   * return `nullopt`.
   */
  std::optional<uint8_t> ms_level() const { return ms_level_; }

  /**
   * Return the delta model need for decoding null marking in the
   * signals file.
   */
  const std::vector<double>& delta_model() const { return delta_model_; }

  /**
   * The scan time in minutes.
   */
  std::optional<double> scan_time() const { return scan_time_; };

private:
  std::optional<uint8_t> ms_level_;
  std::optional<double> scan_time_;
  std::vector<double> delta_model_;
};

} // namespace MzPeak::Metadata
