// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmRTLocal.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/ConsensusMapNormalizerAlgorithmMedian.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/MATH/StatisticFunctions.h>

#include <algorithm>
#include <cmath>
#include <iomanip>

using namespace std;

namespace OpenMS
{
  namespace
  {
    /// One feature handle with a positive intensity.
    struct Cell
    {
      Size feature = 0;        ///< index of the consensus feature
      Size handle = 0;         ///< position of the handle within the consensus feature
      Size map = 0;            ///< map index
      double rt = 0.0;
      double log2_int = 0.0;
      double offset = 0.0;     ///< current log2 offset
      double residual = 0.0;   ///< log2 intensity minus the feature's level (reference cells only)
      bool reference = false;
    };

    /// Median of the residuals of the reference cells sorted[lo, hi), leaving out position @p skip.
    double windowMedian(const vector<Cell>& cells, const vector<Size>& sorted, Size lo, Size hi, Size skip, vector<double>& buffer)
    {
      buffer.clear();
      for (Size i = lo; i < hi; ++i)
      {
        if (i != skip) buffer.push_back(cells[sorted[i]].residual);
      }
      if (buffer.empty()) return 0.0;
      return Math::median(buffer.begin(), buffer.end());
    }
  }

  vector<vector<double>> ConsensusMapNormalizerAlgorithmRTLocal::computeOffsets(const ConsensusMap& map, const Settings& settings,
                                                                                const std::string& acc_filter, const std::string& desc_filter)
  {
    if (settings.neighbors == 0 || settings.rounds == 0 || settings.min_maps == 0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "'neighbors', 'rounds' and 'min_maps' must be at least 1.");
    }

    const Size n_maps = map.getColumnHeaders().size();
    vector<vector<double>> offsets(map.size());
    vector<Cell> cells;
    vector<vector<Size>> cells_of_feature(map.size());

    // collect the cells; a feature with intensities in at least min_maps maps that passes the filters is a reference feature
    for (Size f = 0; f < map.size(); ++f)
    {
      const ConsensusFeature& cf = map[f];
      offsets[f].assign(cf.getFeatures().size(), 0.0);
      Size handle = 0;
      Size n_feature_maps = 0; // handles are sorted by map index: count the distinct ones
      for (const FeatureHandle& fh : cf.getFeatures())
      {
        if (fh.getIntensity() > 0.0)
        {
          if (fh.getMapIndex() >= n_maps)
          {
            throw Exception::ElementNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                             "map index " + std::to_string(fh.getMapIndex()) + " has no column header");
          }
          Cell c;
          c.feature = f;
          c.handle = handle;
          c.map = fh.getMapIndex();
          c.rt = fh.getRT();
          c.log2_int = std::log2(double(fh.getIntensity()));
          if (cells_of_feature[f].empty() || cells[cells_of_feature[f].back()].map != c.map) ++n_feature_maps;
          cells_of_feature[f].push_back(cells.size());
          cells.push_back(c);
        }
        ++handle;
      }
      if (n_feature_maps >= settings.min_maps &&
          ConsensusMapNormalizerAlgorithmMedian::passesFilters_(map.begin() + f, map, acc_filter, desc_filter))
      {
        for (Size c : cells_of_feature[f]) cells[c].reference = true;
      }
    }

    // reference cells of every map, sorted by retention time (ties: in consensus feature order)
    vector<vector<Size>> reference_of_map(n_maps);
    for (Size c = 0; c < cells.size(); ++c)
    {
      if (cells[c].reference) reference_of_map[cells[c].map].push_back(c);
    }
    vector<bool> normalized_map(n_maps, false);
    Size n_normalized_maps = 0;
    for (Size m = 0; m < n_maps; ++m)
    {
      auto& ref = reference_of_map[m];
      std::stable_sort(ref.begin(), ref.end(), [&cells](Size a, Size b) { return cells[a].rt < cells[b].rt; });
      if (ref.size() >= settings.min_reference)
      {
        normalized_map[m] = true;
        ++n_normalized_maps;
      }
    }
    if (n_normalized_maps < 2)
    {
      OPENMS_LOG_WARN << "RT-local normalization: fewer than two maps have at least " << settings.min_reference
                      << " reference features. The maps are not normalized." << endl;
      return offsets;
    }

    // backfitting
    vector<double> buffer;
    buffer.reserve(2 * settings.neighbors);
    for (Size round = 0; round < settings.rounds; ++round)
    {
      for (Size f = 0; f < map.size(); ++f)
      {
        const auto& own = cells_of_feature[f];
        if (own.empty() || !cells[own.front()].reference) continue;
        double level = 0.0;
        for (Size c : own) level += cells[c].log2_int - cells[c].offset;
        level /= double(own.size());
        for (Size c : own) cells[c].residual = cells[c].log2_int - level;
      }
      for (Size m = 0; m < n_maps; ++m)
      {
        if (!normalized_map[m]) continue;
        const auto& ref = reference_of_map[m];
        const Size n = ref.size();
        for (Size i = 0; i < n; ++i)
        {
          const Size lo = (i >= settings.neighbors) ? i - settings.neighbors : 0;
          const Size hi = std::min(n, i + settings.neighbors + 1);
          cells[ref[i]].offset = windowMedian(cells, ref, lo, hi, i, buffer);
        }
      }
    }

    // other cells: median residual of the reference cells of their map around their retention time
    for (Cell& c : cells)
    {
      if (c.reference || !normalized_map[c.map]) continue;
      const auto& ref = reference_of_map[c.map];
      const Size pos = std::lower_bound(ref.begin(), ref.end(), c.rt,
                                        [&cells](Size a, double rt) { return cells[a].rt < rt; }) - ref.begin();
      const Size lo = (pos >= settings.neighbors) ? pos - settings.neighbors : 0;
      const Size hi = std::min(ref.size(), pos + settings.neighbors);
      c.offset = windowMedian(cells, ref, lo, hi, ref.size(), buffer);
    }

    // report and return
    vector<vector<double>> offsets_of_map(n_maps);
    for (const Cell& c : cells)
    {
      offsets[c.feature][c.handle] = c.offset;
      offsets_of_map[c.map].push_back(c.offset);
    }
    for (Size m = 0; m < n_maps; ++m)
    {
      if (!normalized_map[m])
      {
        OPENMS_LOG_INFO << "RT-local normalization: map " << m << " has " << reference_of_map[m].size()
                        << " reference features (fewer than " << settings.min_reference << ") and is not normalized." << endl;
        continue;
      }
      auto& o = offsets_of_map[m];
      const double sd = (o.size() > 1) ? Math::sd(o.begin(), o.end()) : 0.0;
      const double median = Math::median(o.begin(), o.end());
      OPENMS_LOG_INFO << "RT-local normalization: map " << m << ": " << reference_of_map[m].size() << " reference features, median log2 offset "
                      << std::fixed << std::setprecision(3) << median << ", sd along RT " << sd << std::defaultfloat << endl;
    }
    return offsets;
  }

  void ConsensusMapNormalizerAlgorithmRTLocal::normalizeMaps(ConsensusMap& map, const Settings& settings, const std::string& acc_filter, const std::string& desc_filter)
  {
    const vector<vector<double>> offsets = computeOffsets(map, settings, acc_filter, desc_filter);
    for (Size f = 0; f < map.size(); ++f)
    {
      Size handle = 0;
      for (const FeatureHandle& fh : map[f].getFeatures())
      {
        const double o = offsets[f][handle++];
        if (o != 0.0)
        {
          fh.asMutable().setIntensity(fh.getIntensity() * std::exp2(-o));
        }
      }
    }
  }

} // namespace OpenMS
