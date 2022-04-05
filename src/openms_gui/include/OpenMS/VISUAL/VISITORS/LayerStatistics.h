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
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/VISUAL/INTERFACES/IPeptideIds.h>

#include <array>
#include <map>
#include <string>
#include <variant>

namespace OpenMS
{
  class MetaInfoInterface;
  class ConsensusMap;
  class FeatureMap;
  
  /**
     @brief Struct representing the statistics about a set of values

     Min and max are only useful if count > 0
  */
  template <typename VALUE_TYPE>
  struct RangeStats
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

  using RangeStatsInt = RangeStats<int>;
  using RangeStatsDouble = RangeStats<double>;
  using RangeStatsVariant = std::variant<RangeStatsInt, RangeStatsDouble>;

  /// a simple counting struct, for non-numerical occurrences of meta-values
  struct StatsCounter
  {
    size_t counter{0};
  };

  /// Where did a statistic come from? Useful for display to user, and for internal dispatch when user requests a more detailed value distribution
  enum class RangeStatsSource
  {
    CORE,  ///< statistic was obtained from a core data structure of the container, e.g. intensity
    METAINFO, ///< statistic was obtained from MetaInfoInterface of container elements, e.g. "FWHM" for FeatureMaps
    ARRAYINFO,  ///< statistic was obtained from Float/IntegerArrays of the container elements, e.g. "IonMobility" for PeakMap
    SIZE_OF_STATSSOURCE
  };
  
  /// Names corresponding to elements of enum RangeStatsSource
  static const std::array<const char*, (size_t)RangeStatsSource::SIZE_OF_STATSSOURCE> StatsSourceNames = {"core statistics", "meta values", "data arrays"};

  /// Origin and name of a statistic.
  struct RangeStatsType
  {
    RangeStatsSource src;
    std::string name;

    bool operator<(const RangeStatsType& rhs) const
    {
      return std::tie(src, name) < std::tie(rhs.src, rhs.name);
    }

    bool operator==(const RangeStatsType& rhs) const
    {
      return src == rhs.src && name == rhs.name;
    }
  };

  /// collection of Min/Max/Avg statistics from different sources. Note: must be sorted, i.e. do not switch to unordered_map!
  using StatsMap = std::map<RangeStatsType, RangeStatsVariant>;
  /// collection of MetaValues which are not numeric (counts only the number of occurrences per metavalue)
  using StatsCounterMap = std::map<std::string, StatsCounter>;

  /**
      @brief Compute summary statistics (count/min/max/avg) about a container, e.g. intensity, charge, meta values, ...
  */
  class OPENMS_GUI_DLLAPI LayerStatistics
  {
  public:
    /// get all range statistics, any of which can then be plugged into getDistribution()
    const StatsMap& getRangeStatistics() const
    {
      return overview_range_data_;
    }

    /// obtain count statistics for all meta values which are not numerical
    const StatsCounterMap& getCountStatistics() const
    {
      return overview_count_data_;
    }

    /**
       @brief After computing the overview statistic, you can query a concrete distribution by giving the name of the statistic
       @param which Distribution based on which data? 
       @param number_of_bins Number of histogram bins (equally spaced within [min,max] of the distribution)
       @return The distribution
       @throws Exception::InvalidValue if @p which is not a valid overview statistic for the underlying data
    */
    virtual Math::Histogram<> getDistribution(const RangeStatsType& which, const UInt number_of_bins = 500) const = 0;


  protected:
    /// compute the range and count statistics. Call this method in the Ctor of derived classes.
    virtual void computeStatistics_() = 0;
    /// Brings the meta values of one @p meta_interface (a peak or feature) into the statistics
    void bringInMetaStats_(const MetaInfoInterface* meta_interface);

    StatsMap overview_range_data_; ///< data on numerical values computed during getOverviewStatistics 
    StatsCounterMap overview_count_data_; ///< count data on non-numerical values computed during getOverviewStatistics
  };

  /**
    @brief Computes statistics and distributions for a PeakMap    
  */
  class OPENMS_GUI_DLLAPI LayerStatisticsPeakMap
   : public LayerStatistics
  {
  public:
    LayerStatisticsPeakMap(const PeakMap& pm);
 
    Math::Histogram<> getDistribution(const RangeStatsType& which, const UInt number_of_bins) const override;

  private:
    void computeStatistics_() override;
    const PeakMap* pm_; ///< internal reference to a PeakMap -- make sure it does not go out of
                        ///< scope while using this class
  };

  /**
    @brief Computes statistics and distributions for a PeakMap
  */
  class OPENMS_GUI_DLLAPI LayerStatisticsFeatureMap : public LayerStatistics
  {
  public:
    LayerStatisticsFeatureMap(const FeatureMap& fm);

    Math::Histogram<> getDistribution(const RangeStatsType& which,
                                      const UInt number_of_bins) const override;

  private:
    void computeStatistics_() override;
    const FeatureMap* fm_; ///< internal reference to a FeatureMap -- make sure it does not go out of
                           ///< scope while using this class
  };

  /**
    @brief Computes statistics and distributions for a PeakMap
  */
  class OPENMS_GUI_DLLAPI LayerStatisticsConsensusMap : public LayerStatistics
  {
  public:
    LayerStatisticsConsensusMap(const ConsensusMap& cm);

    Math::Histogram<> getDistribution(const RangeStatsType& which,
                                      const UInt number_of_bins) const override;

  private:
    void computeStatistics_() override;
    const ConsensusMap* cm_; ///< internal reference to a PeakMap -- make sure it does not go out of
                             ///< scope while using this class
  };

  /**
    @brief Computes statistics and distributions for a vector<PeptideIdentifications>
  */
  class OPENMS_GUI_DLLAPI LayerStatisticsIdent : public LayerStatistics
  {
  public:
    LayerStatisticsIdent(const IPeptideIds::PepIds& cm);

    Math::Histogram<> getDistribution(const RangeStatsType& which,
                                      const UInt number_of_bins) const override;

  private:
    void computeStatistics_() override;
    const IPeptideIds::PepIds* ids_; ///< internal reference to a PeptideIds -- make sure it does not
                                     ///< go out of scope while using this class
  };
  
} // namespace OpenMS