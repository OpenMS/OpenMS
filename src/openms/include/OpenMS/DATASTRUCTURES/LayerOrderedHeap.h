// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Christie Mathews $
// --------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <functional>
#include <vector>

namespace OpenMS
{

  /**
    @brief Layer-ordered heap for O(n) construction and efficient partial sorting

    Implements the layer-ordered heap (LOH) data structure from:
    "Performing Selection on a Monotonic Function in Lieu of Sorting
    Using Layer-Ordered Heaps" by Lucke, Pennington, Kreitzberg, Kall, Serang
    (J. Proteome Research 2021, DOI: 10.1021/acs.jproteome.0c00711)

    The primary use case is FDR+filter, walking layers to find where FDR crosses a threshold without fully sorting.

    @ingroup Datastructures
  */
  template <typename T, typename Compare = std::less<T>>
  class LayerOrderedHeap
  {
  public:
    //single layer [begin, end) of elements
    struct Layer
    {
      size_t begin;
      size_t end;
      T min_val;
      T max_val;

      size_t size() const { return end - begin; }
    };

    LayerOrderedHeap() = default;

    /**
      @brief Construct LOH from @p data in-place; @p data is rearranged into layer order

      @param[in,out] data Values to organize into layers; rearranged in-place into layer order
      @param[in] comp Comparator defining the ordering (smallest-first by default)
    */
    explicit LayerOrderedHeap(std::vector<T>& data, Compare comp = Compare())
      : n_(data.size()), comp_(comp)
    {
      if (n_ == 0) return;
      buildLayers_(data);
    }

    /// @brief Number of layers
    size_t numLayers() const { return layers_.size(); }

    /// @brief Total number of elements across all layers
    size_t size() const { return n_; }

    /// @brief Whether the heap holds no elements
    bool empty() const { return n_ == 0; }

    /**
      @brief Access the @p i-th layer (ordered from best to worst by @p Compare)

      @param[in] i Layer index in [0, numLayers())
      @return The requested layer's [begin, end) range and min/max metadata
    */
    const Layer& layer(size_t i) const { return layers_[i]; }

    /// @brief All layers, ordered from best to worst by @p Compare
    const std::vector<Layer>& layers() const { return layers_; }

    /**
      @brief Find the layer index whose range contains the @p k-th smallest element

      @param[in] k Zero-based rank in [0, size())
      @return Index of the layer containing rank @p k
    */
    size_t layerForRank(size_t k) const
    {
      assert(k < n_);
      for (size_t i = 0; i < layers_.size(); ++i)
      {
        if (k < layers_[i].end)
          return i;
      }
      return layers_.size() - 1;
    }

  private:
    /**
      @brief Build the layer structure using O(n) partitioning via std::nth_element
      Layers grow geometrically (eg: sizes 1, 2, 4, 8, ...), partition from largest to smallest layer, peeling off the top elements each time
    */
    void buildLayers_(std::vector<T>& data)
    {
      std::vector<size_t> layer_sizes;
      size_t remaining = n_;
      size_t next_size = 1;
      while (remaining > 0)
      {
        size_t actual_size = std::min(next_size, remaining);
        layer_sizes.push_back(actual_size);
        remaining -= actual_size;
        next_size *= 2;
      }

      // build layers from largest to smallest
      size_t active_end = n_;
      std::vector<std::pair<size_t, size_t>> layer_bounds;
      layer_bounds.reserve(layer_sizes.size());

      for (int i = static_cast<int>(layer_sizes.size()) - 1; i >= 0; --i)
      {
        size_t ls = layer_sizes[i];

        if (ls == active_end)
        {
          layer_bounds.push_back({0, ls});
          break;
        }
        size_t split_pos = active_end - ls;

        std::nth_element(data.begin(), data.begin() + split_pos,
                         data.begin() + active_end, comp_);

        layer_bounds.push_back({split_pos, active_end});
        active_end = split_pos;
      }

      std::reverse(layer_bounds.begin(), layer_bounds.end());

      // build layer metadata
      layers_.reserve(layer_bounds.size());
      for (size_t idx = 0; idx < layer_bounds.size(); ++idx)
      {
        size_t b = layer_bounds[idx].first;
        size_t e = layer_bounds[idx].second;
        Layer lay;
        lay.begin = b;
        lay.end = e;
        auto mm = std::minmax_element(data.begin() + b, data.begin() + e, comp_);
        lay.min_val = *mm.first;
        lay.max_val = *mm.second;
        layers_.push_back(lay);
      }
    }

    size_t n_ = 0;
    Compare comp_;
    std::vector<Layer> layers_;
  };

} // namespace OpenMS
