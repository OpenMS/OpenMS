/*

This file is part of the mzpeak project.  It is subject to the license
specified in the LICENSE file which can be found in the top-level
directory of this repository.

*/

#include "mzpeak/chromatogram.h"

#include <algorithm>
#include <ranges>
#include <vector>

#include "mzpeak/exception.h"
#include "mzpeak/schema/psi/array_type.h"
#include "mzpeak/util/delta_estimator.h"

namespace MzPeak {

/******************************************************************************/
Chromatogram::Chromatogram(
    uint64_t index,
    std::shared_ptr<Data::Signals> data,
    const std::vector<Data::ArrayIndex::Dimension>& dims,
    std::unique_ptr<Util::Slice> slice,
    std::shared_ptr<const std::map<uint64_t, ChromatogramMetadata>> md_map)
    : index_(index)
    , decoder_(std::move(data), std::move(slice), Util::DeltaEstimator<double>({}))
    , time_()
    , intensity_()
    , md_map_(std::move(md_map))
{
  // Refuse a role stored in two physical types rather than concatenating them.
  // The array index groups on data type, so a file storing intensity as both
  // float32 and float64 presents two Intensity dimensions; appending both makes
  // one vector of twice the length that looks entirely reasonable.
  auto only_one = [&dims](Schema::PSI::ArrayType role, const char* what) {
    const std::size_t count = std::ranges::count_if(
        dims, [role](const auto& d) { return d.array_type == role; });
    if (count > 1) {
      throw ParquetError(std::string("chromatogram: the ") + what +
                         " array is stored in more than one physical type; "
                         "this reader cannot merge them");
    }
  };
  only_one(Schema::PSI::ArrayType::RelativeTimeOffset, "time");
  only_one(Schema::PSI::ArrayType::Intensity, "intensity");

  for (auto& dim : dims) {
    if (dim.array_type == Schema::PSI::ArrayType::RelativeTimeOffset) {
      decoder_.decimal(dim, time_);

      // Convert to seconds using the unit the file declares, rather than
      // assuming minutes.  Minutes is what every observed file uses and all the
      // specification asks for, so an unrecognised unit means the file is
      // saying something this reader has not been taught -- and silently
      // treating it as minutes would scale every value by 60.
      // An ABSENT unit is refused rather than assumed to be minutes.  The
      // array index schema requires a unit, so absence means either a
      // malformed file or -- reachable in practice -- a coalesced dimension
      // whose columns disagreed, which is exactly when guessing is worst:
      // assuming minutes for a seconds-valued array is a silent 60x error.
      const std::string& unit = decoder_.unit_of(dim);
      if (unit == "UO:0000031") { // minute
        for (auto& t : time_)
          t *= 60.0;
      } else if (unit == "UO:0000010") { // second
        // Already the unit this API reports.
      } else if (unit.empty()) {
        throw ParquetError(
            "chromatogram " + std::to_string(index_) +
            ": the time array declares no unit, or its columns declare "
            "conflicting ones, so it cannot be converted to seconds");
      } else {
        throw ParquetError("chromatogram " + std::to_string(index_) +
                           ": time is in unit '" + unit +
                           "', which this reader cannot convert to seconds");
      }
    } else if (dim.array_type == Schema::PSI::ArrayType::Intensity) {
      decoder_.decimal(dim, intensity_);
      intensity_unit_ = decoder_.unit_of(dim);
    }
  }
}

/******************************************************************************/
const std::vector<Chromatogram::time_type>& Chromatogram::time() const
{
  return time_;
}

/******************************************************************************/
const std::vector<Chromatogram::intensity_type>& Chromatogram::intensity() const
{
  return intensity_;
}

/******************************************************************************/
const std::string& Chromatogram::intensity_unit() const { return intensity_unit_; }

/******************************************************************************/
const ChromatogramMetadata& Chromatogram::metadata() const
{
  static const ChromatogramMetadata empty{};
  if (!md_map_) return empty;
  auto it = md_map_->find(index_);
  return it == md_map_->end() ? empty : it->second;
}

} // namespace MzPeak
