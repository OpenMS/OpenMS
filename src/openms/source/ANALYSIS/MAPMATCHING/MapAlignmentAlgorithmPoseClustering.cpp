// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/ConversionHelper.h>

#include <OpenMS/CONCEPT/Constants.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

namespace OpenMS
{

  MapAlignmentAlgorithmPoseClustering::MapAlignmentAlgorithmPoseClustering() :
    DefaultParamHandler("MapAlignmentAlgorithmPoseClustering"), 
    ProgressLogger(), max_num_peaks_considered_(0)
  {
    defaults_.insert("superimposer:", PoseClusteringAffineSuperimposer().getParameters());
    defaults_.insert("pairfinder:", StablePairFinder().getParameters());
    defaults_.setValue("max_num_peaks_considered", 1000, "The maximal number of peaks/features to be considered per map. To use all, set to '-1'.");
    defaults_.setMinInt("max_num_peaks_considered", -1);

    defaults_.setValue("peak_selection:method", "spread",
                       "How the 'max_num_peaks_considered' MS1 peaks of a peak map (e.g. mzML input) are selected. "
                       "'spread': at most one peak per analyte (the most intense one of its elution profile and isotope pattern), "
                       "spread over the RT and m/z range of the map. "
                       "'intensity': the most intense peaks (the selection of OpenMS 3.6 and earlier); one intense analyte contributes many peaks "
                       "(its elution profile and isotope pattern). Feature maps always use their most intense features.");
    defaults_.setValidStrings("peak_selection:method", {"spread", "intensity"});
    defaults_.setValue("peak_selection:rt_window", 20.0,
                       "('spread' only) A peak within this RT distance (in seconds) of a selected peak at the same m/z, "
                       "or at an isotope position of it, belongs to the same analyte and is not selected. Use about the width of a chromatographic peak.");
    defaults_.setMinFloat("peak_selection:rt_window", 0.0);
    defaults_.setValue("peak_selection:mz_tolerance", 10.0, "('spread' only) m/z tolerance (in ppm) for the same m/z and for isotope positions.");
    defaults_.setMinFloat("peak_selection:mz_tolerance", 0.0);
    defaults_.setValue("peak_selection:max_charge", 4,
                       "('spread' only) Isotope positions m/z +/- k * 1.0033548378 / z (k = 1..3) are suppressed for the charges z = 1..max_charge. 0 suppresses no isotopes.");
    defaults_.setMinInt("peak_selection:max_charge", 0);
    defaults_.setValue("peak_selection:rt_bins", 10, "('spread' only) Number of RT bins over which the selected peaks are spread.");
    defaults_.setMinInt("peak_selection:rt_bins", 1);
    defaults_.setMaxInt("peak_selection:rt_bins", 1000);
    defaults_.setValue("peak_selection:mz_bins", 10, "('spread' only) Number of m/z bins over which the selected peaks are spread.");
    defaults_.setMinInt("peak_selection:mz_bins", 1);
    defaults_.setMaxInt("peak_selection:mz_bins", 1000);
    defaults_.setSectionDescription("peak_selection", "Selection of the MS1 peaks of peak maps");

    defaultsToParam_();
  }

  void MapAlignmentAlgorithmPoseClustering::updateMembers_()
  {
    superimposer_.setParameters(param_.copy("superimposer:", true));
    superimposer_.setLogType(getLogType());

    pairfinder_.setParameters(param_.copy("pairfinder:", true));
    pairfinder_.setLogType(getLogType());

    max_num_peaks_considered_ = param_.getValue("max_num_peaks_considered");

    spread_peaks_ = (param_.getValue("peak_selection:method").toString() == "spread");
    spread_rt_window_ = param_.getValue("peak_selection:rt_window");
    spread_mz_tolerance_ppm_ = param_.getValue("peak_selection:mz_tolerance");
    spread_max_charge_ = static_cast<Size>(int(param_.getValue("peak_selection:max_charge")));
    spread_rt_bins_ = static_cast<Size>(int(param_.getValue("peak_selection:rt_bins")));
    spread_mz_bins_ = static_cast<Size>(int(param_.getValue("peak_selection:mz_bins")));
  }

  MapAlignmentAlgorithmPoseClustering::~MapAlignmentAlgorithmPoseClustering() = default;

  namespace
  {
    // most intense first; equal intensities are ordered by RT, then m/z, so that the selection and its order do not
    // depend on the standard library. A strict weak ordering only without NaN intensities (rejected by the callers).
    template <typename T>
    bool moreIntense(const T& left, const T& right)
    {
      if (left.getIntensity() != right.getIntensity()) return left.getIntensity() > right.getIntensity();
      if (left.getRT() != right.getRT()) return left.getRT() < right.getRT();
      return left.getMZ() < right.getMZ();
    }

    void throwIfNaN(float intensity, const char* what)
    {
      if (std::isnan(intensity))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      std::string(what) + " with a NaN intensity cannot be aligned.", "NaN");
      }
    }

    // NaN or infinite positions break the ordering and the RT x m/z binning
    void throwIfNotFinite(double value, const char* what)
    {
      if (!std::isfinite(value))
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      std::string(what) + " cannot be aligned.", std::to_string(value));
      }
    }
  } // namespace

  void MapAlignmentAlgorithmPoseClustering::selectPeaks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output, Size n)
  {
    output.clear(true);

    // n is capped by the number of MS1 peaks (not by map.getSize(), which also counts the peaks of other MS levels
    // and all chromatogram points)
    Size ms1_peaks = 0;
    for (const MSSpectrum& spectrum : map)
    {
      if (spectrum.getMSLevel() == 1)
      {
        ms1_peaks += spectrum.size();
      }
    }
    n = std::min(n, ms1_peaks);

    // keep the n most intense peaks in a heap whose front is the least intense of them, instead of copying all peaks
    const auto more_intense = [](const Peak2D& left, const Peak2D& right) { return moreIntense(left, right); };
    std::vector<Peak2D> selected;
    selected.reserve(n);
    for (const MSSpectrum& spectrum : map)
    {
      if (spectrum.getMSLevel() != 1)
      {
        continue;
      }
      if (!spectrum.empty())
      {
        throwIfNotFinite(spectrum.getRT(), "An MS1 spectrum with a non-finite RT");
      }
      for (const Peak1D& peak : spectrum)
      {
        throwIfNaN(peak.getIntensity(), "An MS1 peak");
        throwIfNotFinite(peak.getMZ(), "An MS1 peak with a non-finite m/z");
        if (n == 0)
        {
          continue;
        }
        Peak2D candidate;
        candidate.setRT(spectrum.getRT());
        candidate.setMZ(peak.getMZ());
        candidate.setIntensity(peak.getIntensity());
        if (selected.size() < n)
        {
          selected.push_back(candidate);
          std::push_heap(selected.begin(), selected.end(), more_intense);
        }
        else if (more_intense(candidate, selected.front()))
        {
          std::pop_heap(selected.begin(), selected.end(), more_intense);
          selected.back() = candidate;
          std::push_heap(selected.begin(), selected.end(), more_intense);
        }
      }
    }
    std::sort_heap(selected.begin(), selected.end(), more_intense); // most intense first

    output.reserve(selected.size());
    for (Size element_index = 0; element_index < selected.size(); ++element_index)
    {
      output.push_back(ConsensusFeature(map_index, selected[element_index], element_index));
    }

    output.getColumnHeaders()[map_index].size = selected.size();
    output.updateRanges();
  }

  void MapAlignmentAlgorithmPoseClustering::selectSpreadPeaks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output, Size n,
                                                               double rt_window, double mz_tolerance_ppm, Size max_charge,
                                                               Size rt_bins, Size mz_bins)
  {
    output.clear(true);

    // pass 1: RT and m/z range of the MS1 peaks (spans the grid)
    Size ms1_peaks = 0;
    double rt_min = std::numeric_limits<double>::max(), rt_max = std::numeric_limits<double>::lowest();
    double mz_min = std::numeric_limits<double>::max(), mz_max = std::numeric_limits<double>::lowest();
    for (const MSSpectrum& spectrum : map)
    {
      if (spectrum.getMSLevel() != 1 || spectrum.empty())
      {
        continue;
      }
      throwIfNotFinite(spectrum.getRT(), "An MS1 spectrum with a non-finite RT");
      for (const Peak1D& peak : spectrum)
      {
        throwIfNaN(peak.getIntensity(), "An MS1 peak");
        throwIfNotFinite(peak.getMZ(), "An MS1 peak with a non-finite m/z");
        mz_min = std::min(mz_min, peak.getMZ());
        mz_max = std::max(mz_max, peak.getMZ());
      }
      ms1_peaks += spectrum.size();
      rt_min = std::min(rt_min, spectrum.getRT());
      rt_max = std::max(rt_max, spectrum.getRT());
    }
    n = std::min(n, ms1_peaks);
    rt_bins = std::max<Size>(rt_bins, 1);
    mz_bins = std::max<Size>(mz_bins, 1);
    const Size cells = rt_bins * mz_bins;
    // bin of a value in [min, max]; values outside fall into the first or last bin
    const auto binOf = [](double value, double min, double max, Size bins)
    {
      if (max <= min || value <= min) return Size(0);
      const double position = (value - min) / (max - min) * double(bins);
      return position >= double(bins) ? bins - 1 : static_cast<Size>(position);
    };
    const auto cellOf = [&](double rt, double mz)
    {
      return binOf(rt, rt_min, rt_max, rt_bins) * mz_bins + binOf(mz, mz_min, mz_max, mz_bins);
    };
    const auto more_intense = [](const Peak2D& left, const Peak2D& right) { return moreIntense(left, right); };

    // accepted peaks of the same analyte are found through a hash grid with cells at least as large as the RT window
    // and the m/z tolerance
    const double rt_cell = rt_window > 0.0 ? rt_window : 1.0;
    const double mz_cell = std::max(std::abs(mz_max) * mz_tolerance_ppm * 1e-6, 1e-6);
    const auto key = [&](Int64 rt_key, Int64 mz_key) { return (UInt64(rt_key) << 32) ^ UInt64(mz_key); };
    std::vector<double> offsets{0.0};
    constexpr Size max_isotope = 3;
    for (Size z = 1; z <= max_charge; ++z)
    {
      for (Size k = 1; k <= max_isotope; ++k)
      {
        offsets.push_back(double(k) * Constants::C13C12_MASSDIFF_U / double(z));
        offsets.push_back(-double(k) * Constants::C13C12_MASSDIFF_U / double(z));
      }
    }
    std::vector<Peak2D> accepted;
    std::unordered_map<UInt64, std::vector<Size>> grid;
    const auto sameAnalyteAccepted = [&](const Peak2D& candidate)
    {
      const auto rt_low = static_cast<Int64>(std::floor((candidate.getRT() - rt_window) / rt_cell));
      const auto rt_high = static_cast<Int64>(std::floor((candidate.getRT() + rt_window) / rt_cell));
      for (double offset : offsets)
      {
        const double mz = candidate.getMZ() + offset;
        const double tolerance = std::abs(mz) * mz_tolerance_ppm * 1e-6;
        const auto mz_low = static_cast<Int64>(std::floor((mz - tolerance) / mz_cell));
        const auto mz_high = static_cast<Int64>(std::floor((mz + tolerance) / mz_cell));
        for (Int64 rt_key = rt_low; rt_key <= rt_high; ++rt_key)
        {
          for (Int64 mz_key = mz_low; mz_key <= mz_high; ++mz_key)
          {
            const auto it = grid.find(key(rt_key, mz_key));
            if (it == grid.end()) continue;
            for (Size index : it->second)
            {
              if (std::abs(accepted[index].getRT() - candidate.getRT()) <= rt_window &&
                  std::abs(accepted[index].getMZ() - mz) <= tolerance)
              {
                return true;
              }
            }
          }
        }
      }
      return false;
    };
    // Each grid cell keeps its most intense peaks as candidates, at first a multiple of its share of n (one analyte
    // contributes many peaks), so that memory stays proportional to n. A cell can need more (see below); its limit
    // is then raised and the selection repeated.
    constexpr Size candidates_per_selected_peak = 50;
    std::vector<Size> limit(cells, std::min(ms1_peaks, ((n + cells - 1) / cells) * candidates_per_selected_peak));
    std::vector<std::pair<Size, Size>> ranked; // (rank in cell, index into accepted)
    while (n > 0)
    {
      // pass 2: the most intense candidates of each cell, in heaps whose front is the least intense candidate
      std::vector<std::vector<Peak2D>> pools(cells);
      std::vector<bool> truncated(cells, false);
      for (const MSSpectrum& spectrum : map)
      {
        if (spectrum.getMSLevel() != 1)
        {
          continue;
        }
        for (const Peak1D& peak : spectrum)
        {
          Peak2D candidate;
          candidate.setRT(spectrum.getRT());
          candidate.setMZ(peak.getMZ());
          candidate.setIntensity(peak.getIntensity());
          const Size cell = cellOf(candidate.getRT(), candidate.getMZ());
          std::vector<Peak2D>& pool = pools[cell];
          if (pool.size() < limit[cell])
          {
            pool.push_back(candidate);
            std::push_heap(pool.begin(), pool.end(), more_intense);
            continue;
          }
          truncated[cell] = true;
          if (more_intense(candidate, pool.front()))
          {
            std::pop_heap(pool.begin(), pool.end(), more_intense);
            pool.back() = candidate;
            std::push_heap(pool.begin(), pool.end(), more_intense);
          }
        }
      }
      std::vector<Peak2D> least_intense(cells); // the least intense candidate of each truncated cell
      for (Size cell = 0; cell < cells; ++cell)
      {
        if (truncated[cell])
        {
          least_intense[cell] = pools[cell].front();
        }
      }
      std::vector<Peak2D> candidates;
      for (std::vector<Peak2D>& pool : pools)
      {
        candidates.insert(candidates.end(), pool.begin(), pool.end());
        std::vector<Peak2D>().swap(pool);
      }
      std::sort(candidates.begin(), candidates.end(), more_intense);

      // pass 3: from the most intense candidate down, accept a candidate unless an accepted peak of the same analyte
      // (same m/z or an isotope position, within rt_window) exists
      accepted.clear();
      grid.clear();
      for (const Peak2D& candidate : candidates)
      {
        if (sameAnalyteAccepted(candidate))
        {
          continue;
        }
        grid[key(static_cast<Int64>(std::floor(candidate.getRT() / rt_cell)),
                 static_cast<Int64>(std::floor(candidate.getMZ() / mz_cell)))].push_back(accepted.size());
        accepted.push_back(candidate);
      }

      // round-robin over the grid cells: rank the accepted peaks within their cell (they are in intensity order) and
      // take the n peaks of lowest rank, the more intense first within a rank
      std::vector<Size> accepted_in_cell(cells, 0);
      ranked.clear();
      ranked.reserve(accepted.size());
      for (Size i = 0; i < accepted.size(); ++i)
      {
        ranked.emplace_back(accepted_in_cell[cellOf(accepted[i].getRT(), accepted[i].getMZ())]++, i);
      }
      std::sort(ranked.begin(), ranked.end());
      const bool all_taken = ranked.size() < n;
      ranked.resize(std::min(n, ranked.size()));

      // A truncated cell needs more candidates if
      // (a) it may hold further analytes that the round-robin would take (it reaches the rank of the cell's next
      //     accepted peak): wanted[cell] is the number of accepted peaks the cell should yield, or
      // (b) one of its cut-off peaks may belong to the analyte of a selected peak of another cell: the selected peak
      //     lies within reach (as in pass 3) and is not more intense than the cell's least intense candidate. It is
      //     then e.g. the tail of an analyte whose apex was cut off.
      const Size last_rank = ranked.empty() ? 0 : ranked.back().first;
      std::vector<Size> wanted(cells, 0);
      Size short_cells = 0;
      for (Size cell = 0; cell < cells; ++cell)
      {
        if (truncated[cell] && (all_taken || accepted_in_cell[cell] <= last_rank))
        {
          ++short_cells;
        }
      }
      for (Size cell = 0; cell < cells && short_cells > 0; ++cell)
      {
        if (truncated[cell] && (all_taken || accepted_in_cell[cell] <= last_rank))
        {
          // if all accepted peaks are selected, the missing ones are shared among these cells; otherwise the cell
          // should reach the last selected rank
          wanted[cell] = all_taken ? accepted_in_cell[cell] + (n - ranked.size() + short_cells - 1) / short_cells : last_rank + 1;
        }
      }
      std::vector<bool> cut_suppressor(cells, false);
      const double isotope_span = max_charge > 0 ? double(max_isotope) * Constants::C13C12_MASSDIFF_U : 0.0;
      for (const auto& entry : ranked)
      {
        const Peak2D& peak = accepted[entry.second];
        const Size own_cell = cellOf(peak.getRT(), peak.getMZ());
        const double tolerance = std::abs(peak.getMZ() + isotope_span) * mz_tolerance_ppm * 1e-6;
        const Size rt_first = binOf(peak.getRT() - rt_window, rt_min, rt_max, rt_bins);
        const Size rt_last = binOf(peak.getRT() + rt_window, rt_min, rt_max, rt_bins);
        const Size mz_first = binOf(peak.getMZ() - isotope_span - tolerance, mz_min, mz_max, mz_bins);
        const Size mz_last = binOf(peak.getMZ() + isotope_span + tolerance, mz_min, mz_max, mz_bins);
        for (Size rt_bin = rt_first; rt_bin <= rt_last; ++rt_bin)
        {
          for (Size mz_bin = mz_first; mz_bin <= mz_last; ++mz_bin)
          {
            const Size cell = rt_bin * mz_bins + mz_bin;
            if (cell != own_cell && truncated[cell] && !more_intense(peak, least_intense[cell]))
            {
              cut_suppressor[cell] = true;
            }
          }
        }
      }

      bool repeat = false;
      for (Size cell = 0; cell < cells; ++cell)
      {
        double factor = 0.0;
        if (wanted[cell] > 0)
        {
          // in proportion to the shortfall, with a margin
          factor = std::clamp(1.5 * double(wanted[cell]) / double(std::max<Size>(accepted_in_cell[cell], 1)), 2.0, 16.0);
        }
        else if (cut_suppressor[cell])
        {
          factor = 4.0;
        }
        if (factor > 0.0)
        {
          limit[cell] = std::min(ms1_peaks, static_cast<Size>(std::ceil(double(limit[cell]) * factor)));
          repeat = true;
        }
      }
      if (!repeat)
      {
        break;
      }
    }

    std::vector<Peak2D> selected;
    selected.reserve(ranked.size());
    for (const auto& entry : ranked)
    {
      selected.push_back(accepted[entry.second]);
    }
    std::sort(selected.begin(), selected.end(), more_intense);

    output.reserve(selected.size());
    for (Size element_index = 0; element_index < selected.size(); ++element_index)
    {
      output.push_back(ConsensusFeature(map_index, selected[element_index], element_index));
    }
    output.getColumnHeaders()[map_index].size = selected.size();
    output.updateRanges();
  }

  void MapAlignmentAlgorithmPoseClustering::selectPeakLandmarks_(UInt64 map_index, const PeakMap& map, ConsensusMap& output) const
  {
    if (spread_peaks_)
    {
      selectSpreadPeaks_(map_index, map, output, max_num_peaks_considered_, spread_rt_window_, spread_mz_tolerance_ppm_,
                         spread_max_charge_, spread_rt_bins_, spread_mz_bins_);
    }
    else
    {
      selectPeaks_(map_index, map, output, max_num_peaks_considered_);
    }
  }

  void MapAlignmentAlgorithmPoseClustering::selectFeatures_(UInt64 map_index, const FeatureMap& map, ConsensusMap& output, Size n)
  {
    for (const Feature& feature : map)
    {
      throwIfNaN(feature.getIntensity(), "A feature");
      throwIfNotFinite(feature.getRT(), "A feature with a non-finite RT");
      throwIfNotFinite(feature.getMZ(), "A feature with a non-finite m/z");
    }
    n = std::min(n, map.size());

    // indices of the n most intense features (ties by RT, m/z and index), in map order
    std::vector<Size> selected(map.size());
    std::iota(selected.begin(), selected.end(), 0);
    if (n < map.size())
    {
      std::nth_element(selected.begin(), selected.begin() + n, selected.end(), [&map](Size left, Size right)
      {
        if (moreIntense(map[left], map[right])) return true;
        if (moreIntense(map[right], map[left])) return false;
        return left < right;
      });
      selected.resize(n);
      std::sort(selected.begin(), selected.end());
    }

    // the identifications of the selected features stay, so that the pair finder can use them
    MapConversion::convert(map_index, map, output, selected);
  }

  void MapAlignmentAlgorithmPoseClustering::setReference(const FeatureMap& map)
  {
    selectFeatures_(0, map, reference_, max_num_peaks_considered_);
  }

  void MapAlignmentAlgorithmPoseClustering::setReference(const PeakMap& map)
  {
    selectPeakLandmarks_(0, map, reference_);
  }

  void MapAlignmentAlgorithmPoseClustering::align(const FeatureMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    selectFeatures_(1, map, map_scene, max_num_peaks_considered_);
    align(map_scene, trafo);
  }

  void MapAlignmentAlgorithmPoseClustering::align(const PeakMap& map, TransformationDescription& trafo)
  {
    ConsensusMap map_scene;
    selectPeakLandmarks_(1, map, map_scene);
    align(map_scene, trafo);
  }

  void MapAlignmentAlgorithmPoseClustering::align(const ConsensusMap& map, TransformationDescription& trafo)
  {
    // TODO: move this to updateMembers_? (if ConsensusMap prevails)
    // TODO: why does superimposer work on consensus map???
    const ConsensusMap & map_model = reference_;
    ConsensusMap map_scene = map;

    // run superimposer to find the global transformation
    TransformationDescription si_trafo;
    superimposer_.run(map_model, map_scene, si_trafo);

    // apply transformation to consensus features and contained feature
    // handles
    for (Size j = 0; j < map_scene.size(); ++j)
    {
      // Calculate new RT
      double rt = map_scene[j].getRT();
      rt = si_trafo.apply(rt);
      // Set RT of consensus feature centroid
      map_scene[j].setRT(rt);
      // Set RT of consensus feature handles
      map_scene[j].begin()->asMutable().setRT(rt);
    }

    // run pairfinder to find pairs
    ConsensusMap result;
    // TODO: add another 2-map interface to pairfinder?
    std::vector<ConsensusMap> input(2);
    input[0] = map_model;
    input[1] = map_scene;
    pairfinder_.run(input, result);

    // calculate the local transformation
    si_trafo.invert(); // to undo the transformation applied above
    TransformationDescription::DataPoints data;
    for (ConsensusFeature& cfeature : result)
    {
      if (cfeature.size() == 2) // two matching features
      {
        ConsensusFeature::iterator feat_it = cfeature.begin();
        double y = feat_it->getRT();
        double x = si_trafo.apply((++feat_it)->getRT());
        // one feature should be from the reference map:
        if (feat_it->getMapIndex() != 0)
        {
          data.push_back(make_pair(x, y));
        }
        else
        {
          data.push_back(make_pair(y, x));
        }
      }
    }
    trafo = TransformationDescription(data);
    trafo.fitModel("linear");
  }

} // namespace
