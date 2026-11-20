/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#pragma once

#include <memory>
#include <vector>

namespace MzPeak::Util {
class Parquet;
}

namespace MzPeak::Metadata {

/**
 * The method of precursor-ion selection and activation.
 */
class Scans final {
public:
  /// Data for all scans.
  struct Raw {
    std::vector<double> scan_start_time;
  };

  /// Default constructor;
  Scans();

  /// Constructor.
  Scans(std::unique_ptr<Util::Parquet>, uint64_t);

  /**
   * Return the raw table data for all scans.
   */
  const Raw& raw() const { return raw_; }

private:
  Raw raw_;
};

} // namespace MzPeak::Metadata
