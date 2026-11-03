// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
#include "OpenDIARegionWorkflow.h"

#include "OpenDIARawAccess.h"
#include "OpenDIASparseEvidence.h"

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessQuadMZTransforming.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

namespace OpenDIA
{
namespace
{
  constexpr OpenMS::Size MAX_SPARSE_EVIDENCE_ROWS = 131072;
  constexpr OpenMS::Size MAX_SPARSE_QUERIES = 131072;
  using Coordinate = OpenMS::ChromatogramExtractor::ExtractionCoordinates;
  [[noreturn]] void invalid(const char* message)
  {
    throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, message);
  }

  double integrate(const RawRegion::Scan& scan, const Coordinate& coordinate, const OpenMS::ChromExtractParams& cp, OpenMS::Size& center)
  {
    const auto mz = scan.mz();
    const auto intensity = scan.intensity();
    const auto im = scan.mobility();
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
  // Only native region columns share the raw budget. This is not a process RSS cap.
  std::vector<double> rts(count);
  for (OpenMS::Size i = 0; i < count; ++i)
  {
    rts[i] = raw->getSpectrumMetaById(static_cast<int>(i)).RT;
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
  // The cache keys each region by source, namespace and exact calibration coefficients.
  const RecordContext context {0, 0, 0, calibration ? 1u : 0u};
  std::optional<SparseEvidenceExtractor> sparse;
  bool sparse_bounds_valid = true;
  if (cp.im_extraction_window > 0)
  {
    const double half_im_width = cp.im_extraction_window / 2.0;
    for (const auto& coordinate : coordinates)
    {
      if (coordinate.ion_mobility >= 0 &&
          (!std::isfinite(coordinate.ion_mobility - half_im_width) ||
           !std::isfinite(coordinate.ion_mobility + half_im_width)))
      {
        sparse_bounds_valid = false;
        break;
      }
    }
  }
  if (coordinates.size() <= MAX_SPARSE_QUERIES && sparse_bounds_valid)
  {
    std::vector<ExtractionQuery> selectors;
    std::vector<FragmentQuery> queries;
    selectors.reserve(coordinates.size());
    queries.reserve(coordinates.size());
    for (OpenMS::Size k = 0; k < coordinates.size(); ++k)
    {
      const auto& coordinate = coordinates[k];
      const double half_width = cp.ppm ? coordinate.mz * cp.mz_extraction_window / 2.0 * 1.0e-6 : cp.mz_extraction_window / 2.0;
      const double left = coordinate.mz - half_width;
      const double right = coordinate.mz + half_width;
      ExtractionQuerySpec selector;
      selector.context = context;
      selector.ms_level = 0;
      if (coordinate.rt_end - coordinate.rt_start > 0)
      {
        selector.rt_seconds = QueryInterval {coordinate.rt_start, coordinate.rt_end, EdgePolicy::CLOSED};
      }
      if (cp.im_extraction_window > 0 && coordinate.ion_mobility >= 0)
      {
        const double half_im_width = cp.im_extraction_window / 2.0;
        // RegionWorkflow coordinates use the raw point-mobility numeric scale.
        selector.mobility = MobilityInterval {{coordinate.ion_mobility - half_im_width,
                                                coordinate.ion_mobility + half_im_width,
                                                EdgePolicy::OPEN}, MobilityUnit::INVERSE_REDUCED};
      }
      selector.mz = {left, right, EdgePolicy::OPEN};
      selector.tolerance = {cp.mz_extraction_window, cp.ppm ? ToleranceUnit::PPM : ToleranceUnit::TH};
      selectors.emplace_back(std::move(selector));
      queries.push_back({0, 0, static_cast<RecordId>(k), k, {left, right, EdgePolicy::OPEN}});
    }
    sparse.emplace(context, RecordLimits {coordinates.size(), MAX_SPARSE_EVIDENCE_ROWS},
                   std::move(selectors), std::move(queries));
  }
  std::vector<SparseEvidenceRef> evidence;
  for (OpenMS::Size i = 0; i < count; ++i)
  {
    const auto region = raw->acquireRegion(static_cast<int>(i), context, calibration);
    const auto scan = region.scan(0);
    if (scan.mz().empty()) { continue; }
    if (cp.im_extraction_window > 0 && ! scan.acquisition().point_mobility) { invalid("Requested point mobility is absent"); }
    if (!sparse || !sparse->tryExtract(scan, evidence))
    {
      OpenMS::Size center = 0;
      for (OpenMS::Size k = 0; k < coordinates.size(); ++k)
      {
        const auto& coordinate = coordinates[k];
        if (coordinate.rt_end - coordinate.rt_start > 0 &&
            (rts[i] < coordinate.rt_start || rts[i] > coordinate.rt_end)) { continue; }
        output[k]->getTimeArray()->data.push_back(rts[i]);
        output[k]->getIntensityArray()->data.push_back(integrate(scan, coordinate, cp, center));
      }
      continue;
    }
    std::sort(evidence.begin(), evidence.end(), [](const auto& lhs, const auto& rhs)
    {
      if (lhs.fragment_id != rhs.fragment_id) { return lhs.fragment_id < rhs.fragment_id; }
      return lhs.raw_point_index < rhs.raw_point_index;
    });
    const auto mz = scan.mz();
    const auto intensity = scan.intensity();
    OpenMS::Size center = 0;
    for (OpenMS::Size k = 0; k < coordinates.size(); ++k)
    {
      const auto& coordinate = coordinates[k];
      if (coordinate.rt_end - coordinate.rt_start > 0 &&
          (rts[i] < coordinate.rt_start || rts[i] > coordinate.rt_end)) { continue; }
      output[k]->getTimeArray()->data.push_back(rts[i]);
      const auto query_id = static_cast<RecordId>(k);
      const auto first = std::lower_bound(evidence.begin(), evidence.end(), query_id,
        [](const SparseEvidenceRef& row, RecordId id) { return row.fragment_id < id; });
      const auto last = std::upper_bound(first, evidence.end(), query_id,
        [](RecordId id, const SparseEvidenceRef& row) { return id < row.fragment_id; });
      const double half_width = cp.ppm ? coordinate.mz * cp.mz_extraction_window / 2.0 * 1.0e-6 : cp.mz_extraction_window / 2.0;
      const double left = coordinate.mz - half_width;
      const double right = coordinate.mz + half_width;
      center = static_cast<OpenMS::Size>(std::lower_bound(mz.begin() + center, mz.end(), coordinate.mz) - mz.begin());
      double total = 0;
      const auto center_row = std::lower_bound(first, last, center,
        [](const SparseEvidenceRef& row, OpenMS::Size point_index) { return row.raw_point_index < point_index; });
      if (center < mz.size() && mz[center] > left && mz[center] < right && center_row != last &&
          center_row->raw_point_index == center)
      {
        total += intensity[center];
      }
      for (auto row = center_row; row != first;)
      {
        --row;
        if (mz[row->raw_point_index] > left) { total += intensity[row->raw_point_index]; }
      }
      for (auto row = center_row; row != last; ++row)
      {
        if (row->raw_point_index > center && mz[row->raw_point_index] < right)
        {
          total += intensity[row->raw_point_index];
        }
      }
      output[k]->getIntensityArray()->data.push_back(total);
    }
  }
}
} // namespace OpenDIA
