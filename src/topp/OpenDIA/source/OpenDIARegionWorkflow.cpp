// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include "OpenDIARegionWorkflow.h"

#include "OpenDIARawAccess.h"

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <span>
#include <utility>

namespace OpenDIA
{
namespace
{
  using Coordinate = OpenMS::ChromatogramExtractor::ExtractionCoordinates;
  [[noreturn]] void invalid(const char* message)
  {
    throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
  }

  // Private to the audited read-only scoring path. Unlike RawAccess's detached
  // mutable API, this access shares the same admitted columns used for XICs.
  class ScoringAccess final : public OpenMS::SpectrumAccessTransforming
  {
  public:
    ScoringAccess(std::shared_ptr<RawAccess> raw, std::optional<MZCalibration> calibration) :
      SpectrumAccessTransforming(raw), raw_(std::move(raw)), calibration_(calibration)
    {
    }

    OpenSwath::SpectrumPtr getSpectrumById(int id) override
    {
      return raw_->spectrumForScoring(id, calibration_);
    }

    std::shared_ptr<OpenSwath::ISpectrumAccess> lightClone() const override
    {
      return std::make_shared<ScoringAccess>(std::dynamic_pointer_cast<RawAccess>(raw_->lightClone()), calibration_);
    }

  private:
    std::shared_ptr<RawAccess> raw_;
    std::optional<MZCalibration> calibration_;
  };

  template<typename Intensity>
  double integrate(std::span<const double> mz, std::span<const Intensity> intensity, std::span<const Intensity> im,
                   const Coordinate& coordinate, const OpenMS::ChromExtractParams& cp, OpenMS::Size& center)
  {
    const double half_width = cp.ppm ? coordinate.mz * cp.mz_extraction_window / 2.0 * 1.0e-6 : cp.mz_extraction_window / 2.0;
    const double left = coordinate.mz - half_width, right = coordinate.mz + half_width;
    const bool use_im = cp.im_extraction_window > 0 && coordinate.ion_mobility >= 0;
    const double left_im = coordinate.ion_mobility - cp.im_extraction_window / 2.0;
    const double right_im = coordinate.ion_mobility + cp.im_extraction_window / 2.0;
    center = static_cast<OpenMS::Size>(std::lower_bound(mz.begin() + center, mz.end(), coordinate.mz) - mz.begin());
    double total = 0;
    const auto add = [&](OpenMS::Size p) {
      if (! use_im || (im[p] > left_im && im[p] < right_im)) { total += intensity[p]; }
    };
    // Preserve the reference's center, descending-left, ascending-right addition order.
    // Changing that order can alter downstream peak selection through rounding.
    if (center < mz.size() && mz[center] > left && mz[center] < right) { add(center); }
    for (auto p = center; p > 0;)
    {
      --p;
      if (! (mz[p] > left && mz[p] < right)) { break; }
      add(p);
    }
    for (auto p = center < mz.size() ? center + 1 : center; p < mz.size() && mz[p] > left && mz[p] < right; ++p)
    {
      add(p);
    }
    return total;
  }
} // namespace

OpenSwath::SpectrumAccessPtr RegionWorkflow::spectrumAccessForScoring_(const OpenSwath::SpectrumAccessPtr& input) const
{
  auto source = input;
  std::optional<MZCalibration> calibration;
  if (const auto transformed = std::dynamic_pointer_cast<OpenMS::SpectrumAccessQuadMZTransforming>(source))
  {
    const auto [a, b, c, ppm] = transformed->getTransformation();
    calibration = MZCalibration {a, b, c, ppm};
    source = transformed->getInput();
  }
  if (auto raw = std::dynamic_pointer_cast<RawAccess>(source))
  {
    return std::make_shared<ScoringAccess>(std::move(raw), calibration);
  }
  // Unknown readers or stacked transformations retain their exact behavior.
  return input;
}

void RegionWorkflow::extractChromatograms_(const OpenSwath::SpectrumAccessPtr& input,
                                           std::vector<OpenSwath::ChromatogramPtr>& output,
                                           const std::vector<Coordinate>& coordinates,
                                           const OpenMS::ChromExtractParams& cp)
{
  auto source = input;
  std::optional<MZCalibration> calibration;
  if (const auto transformed = std::dynamic_pointer_cast<OpenMS::SpectrumAccessQuadMZTransforming>(source))
  {
    const auto [a, b, c, ppm] = transformed->getTransformation();
    calibration = MZCalibration {a, b, c, ppm};
    source = transformed->getInput();
  }
  const auto raw = std::dynamic_pointer_cast<RawAccess>(source);
  if (! raw || cp.extraction_function != "tophat")
  {
    OpenMS::OpenSwathWorkflowBase::extractChromatograms_(input, output, coordinates, cp);
    return;
  }
  if (output.size() != coordinates.size()) { invalid("Output and coordinate counts differ"); }
  if (! std::isfinite(cp.mz_extraction_window) || cp.mz_extraction_window < 0 || ! std::isfinite(cp.im_extraction_window))
  {
    invalid("Invalid extraction width");
  }
  for (OpenMS::Size k = 0; k < coordinates.size(); ++k)
  {
    const auto& c = coordinates[k];
    if (! output[k] || ! output[k]->getTimeArray() || ! output[k]->getIntensityArray()) { invalid("Missing output arrays"); }
    if (! std::isfinite(c.mz) || c.mz < 0 || ! std::isfinite(c.rt_start) || ! std::isfinite(c.rt_end) || ! std::isfinite(c.ion_mobility)
        || (k && coordinates[k - 1].mz > c.mz))
    {
      invalid("Invalid or unsorted extraction coordinates");
    }
    const double half_width = cp.ppm ? c.mz * cp.mz_extraction_window / 2.0 * 1.0e-6 : cp.mz_extraction_window / 2.0;
    if (! std::isfinite(c.mz - half_width) || ! std::isfinite(c.mz + half_width)) { invalid("Extraction bounds overflow"); }
  }
  const auto count = raw->getNrSpectra();
  if (count > static_cast<OpenMS::Size>(std::numeric_limits<int>::max())) { invalid("Too many acquisitions"); }
  // Output arrays and temporary RT metadata retain the workflow's existing allocation policy.
  // Only cached scan columns share the raw budget. This is not a process RSS cap.
  std::vector<double> rts(count);
  for (OpenMS::Size i = 0; i < count; ++i)
  {
    rts[i] = raw->getSpectrumMetaById(static_cast<int>(i)).RT;
    if (! std::isfinite(rts[i])) { invalid("Invalid acquisition RT"); }
  }
  auto sorted_rts = rts;
  std::sort(sorted_rts.begin(), sorted_rts.end());
  for (OpenMS::Size k = 0; k < coordinates.size(); ++k)
  {
    const auto& c = coordinates[k];
    const auto points = c.rt_end - c.rt_start > 0 ? std::upper_bound(sorted_rts.begin(), sorted_rts.end(), c.rt_end)
                                                      - std::lower_bound(sorted_rts.begin(), sorted_rts.end(), c.rt_start)
                                                  : static_cast<std::ptrdiff_t>(count);
    auto& times = output[k]->getTimeArray()->data;
    auto& intensities = output[k]->getIntensityArray()->data;
    times.reserve(times.size() + points);
    intensities.reserve(intensities.size() + points);
  }
  const OpenSwath::SpectrumAccessScope scope(OpenSwath::SpectrumAccessScope::Purpose::EXTRACTION);
  // Extraction and scoring use one cache key, including the exact calibration.
  // Once their last owner is gone, the LRU may reclaim these admitted columns.
  for (OpenMS::Size i = 0; i < count; ++i)
  {
    auto scan = raw->acquire(static_cast<int>(i), calibration);
    const bool sorted = scan->validateForExtraction();
    const auto append = [&]<typename Intensity>(std::span<const double> mz, std::span<const Intensity> intensity,
                                               std::span<const Intensity> im)
    {
      if (mz.empty()) { return; }
      if (cp.im_extraction_window > 0 && im.empty()) { invalid("Requested point mobility is absent"); }
      OpenMS::Size center = 0;
      for (OpenMS::Size k = 0; k < coordinates.size(); ++k)
      {
        const auto& c = coordinates[k];
        if (c.rt_end - c.rt_start > 0 && (rts[i] < c.rt_start || rts[i] > c.rt_end)) { continue; }
        output[k]->getTimeArray()->data.push_back(rts[i]);
        output[k]->getIntensityArray()->data.push_back(integrate(mz, intensity, im, c, cp, center));
      }
    };
    if (sorted)
    {
      const auto mobility = scan->mobility();
      append(std::span<const double>(scan->array(0).data), std::span<const double>(scan->array(1).data),
             mobility ? std::span<const double>(mobility->data) : std::span<const double>());
    }
    else
    {
      // Rare non-monotonic transformations retain extraction's stable ordering.
      // Scoring still sees the original transformed order, as in the legacy API.
      scan.reset();
      const auto region = raw->acquireRegion(static_cast<int>(i), {}, calibration).scan(0);
      append(region.mz(), region.intensity(), region.mobility());
    }
  }
}
} // namespace OpenDIA
