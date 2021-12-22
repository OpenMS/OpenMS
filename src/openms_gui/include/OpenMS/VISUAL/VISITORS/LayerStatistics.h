// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <map>
#include <string>
#include <variant>

namespace OpenMS
{
  class MetaInfoInterface;
  class ConsensusMap;
  class FeatureMap;
  
  /**
     @brief Struct representing the statistics about one meta information.

     Min and max are only useful if count > 0
  */
  template <typename VALUE_TYPE>
  struct StatsSummary
  {
    public:  
      void addDataPoint(VALUE_TYPE v)
      {
        ++count_;
        sum_ += v;
        min_ = std::min(min_, v);
        max_ = std::max(max_, v);
      }

      VALUE_TYPE getMin() const
      {
        return min_;
      }

      VALUE_TYPE getMax() const
      {
        return max_;
      }

      size_t getCount() const
      {
        return count_;
      }

      /// get the average value from all calls to addDataPoint()
      double getAvg() const
      {
        return count_ == 0 ? 0 : double(sum_) / count_;
      }

    protected:
      size_t count_{0};
      VALUE_TYPE min_{std::numeric_limits<VALUE_TYPE>::max()}; // init with very high value 
      VALUE_TYPE max_{std::numeric_limits<VALUE_TYPE>::lowest()}; // init with lowest (=negative) value possible
      VALUE_TYPE sum_{0};
  };

  /// a simple counting struct, for non-numerical occurrences of meta-values
  struct StatsCounter
  {
    size_t counter{0};
  };

  using SSInt = StatsSummary<int>;
  using SSDouble = StatsSummary<double>;
  using StatsSummaryVariant = std::variant<SSInt, SSDouble, StatsCounter>;

  /**
      @brief Compute summary statistics (count/min/max/avg) about a layer, e.g. intensity, charge, meta values, ...
  */
  class OPENMS_GUI_DLLAPI LayerStatistics
  {
  public:
    /// Computes the statistics of a peak layer
    void computePeakMapStats(const PeakMap& pm);
    /// Computes the statistics of a feature layer
    void computeFeatureMapStats(const FeatureMap& fm);
    /// Computes the statistics of a consensus feature layer
    void computeConsensusMapStats(const ConsensusMap& cm);

    using StatsMap = std::map<std::string, StatsSummaryVariant>;

    /// retrieve the core statistics of the data structure (usually intensity, maybe charge, quality, etc)
    const StatsMap& getCoreStats() const;

    /// retrieve the statistics of all the meta values, e.g. "FWHM" for FeatureMaps
    /// The result is computed by looking up a MetaInfoInterface. Try not to query repeatedly if performance matters.
    StatsMap getMetaStats() const;

    /// retrieve the statistics of all float and integer data arrays of spectra (only populated for PeakMaps)
    const StatsMap& getArrayStats() const;

  protected:
    /// Brings the meta values of one @p meta_interface (a peak or feature) into the statistics
    void bringInMetaStats_(const MetaInfoInterface* meta_interface);
    
    /// core statistics of the data structure
    StatsMap core_stats_;
    /// Map containing the statistics about all meta information of the peaks/features in the layer
    /// It's key is UInt for performance reasons
    std::map<UInt, StatsSummaryVariant> meta_stats_;
    /// Map containing the statistics about the FloatDataArrays of all spectra in this layer
    StatsMap meta_array_stats_;
  };

} // namespace OpenMS