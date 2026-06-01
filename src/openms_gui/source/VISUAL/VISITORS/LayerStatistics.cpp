// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/MetaInfo.h>

using namespace std;

namespace OpenMS
{
  // helper for visitor pattern with std::visit
  template<class... Ts>
  struct overload : Ts... {
    using Ts::operator()...;
  };
  template<class... Ts>
  overload(Ts...) -> overload<Ts...>;

  struct MinMax {
    double min;
    double max;
  };

  /// extract min,max from a statistic (or throw Exception if stats is not present in @p overview_data)
  MinMax getMinMax(const StatsMap& overview_data, const RangeStatsType& which, const std::string& error_message_container)
  {
    auto overview_stat = overview_data.find(which);
    if (overview_stat == overview_data.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Statistic is not valid for this " + error_message_container,
                                    which.name);
    }
    // getMin/Max from variant
    MinMax result = std::visit(overload {[](const auto& stats) -> MinMax {
                                 return {(double)stats.getMin(), (double)stats.getMax()};
                               }},
        overview_stat->second);
    
    return result;
  }

  /// Computes the statistics of all meta data contained in the FloatDataArray or IntegerDataArray of an MSSpectrum
  template<typename SSType, typename DataArrayType>
  void computeMetaDataArrayStats(const DataArrayType& arrays, StatsMap& stats)
  {
    using VarType = RangeStats<SSType>;

    for (const auto& mda : arrays)
    {
      const RangeStatsType mda_name = {RangeStatsSource::ARRAYINFO, mda.getName()};
      auto it_var = stats.find(mda_name);
      if (it_var == stats.end()) // missing, create it
      {
        it_var = stats.emplace(mda_name, VarType()).first;
      }
      RangeStatsVariant& meta_stats_value = it_var->second;
      auto& target = std::get<VarType>(meta_stats_value);
      for (const auto& value : mda)
      {
        target.addDataPoint(value);
      }
    }
  }

  /// Update the histogram for data of a certain FloatDataArray or IntegerDataArray
  /// of an MSSpectrum
  template<typename DataArrayType>
  void updateHistFromDataArray(const DataArrayType& arrays, const std::string& name, Math::Histogram<>& hist)
  {
    for (const auto& mda : arrays)
    {
      if (name != mda.getName()) continue;
      for (const auto& value : mda)
      {
        hist.inc(value);
      }
    }
  }

  LayerStatisticsPeakMap::LayerStatisticsPeakMap(const PeakMap& pm) 
    : pm_(&pm)
  {
    computeStatistics_();
  }

  void LayerStatisticsPeakMap::computeStatistics_()
  {
    RangeStats<double> stat_intensity;
    for (const auto& spec : *pm_)
    {
      for (const auto& peak : spec)
      {
        stat_intensity.addDataPoint(peak.getIntensity());
      }
      // collect stats about the meta data arrays of this spectrum
      computeMetaDataArrayStats<double>(spec.getFloatDataArrays(), overview_range_data_);
      computeMetaDataArrayStats<int>(spec.getIntegerDataArrays(), overview_range_data_);
    }
    overview_range_data_.emplace(RangeStatsType {RangeStatsSource::CORE, "intensity"}, stat_intensity);
  }

  Math::Histogram<> LayerStatisticsPeakMap::getDistribution(const RangeStatsType& which, const UInt number_of_bins) const
  {
    auto mm = getMinMax(overview_range_data_, which, "PeakMap"); // may throw if unknown statistic
    Math::Histogram<> result(mm.min, mm.max, (mm.max - mm.min) / number_of_bins);
    
    if (which == RangeStatsType{ RangeStatsSource::CORE, "intensity" })
    {
      for (const auto& spec : *pm_)
      {
        for (const auto& peak : spec)
        {
          result.inc(peak.getIntensity());
        }
      }
    }
    else if (which.src == RangeStatsSource::ARRAYINFO)
    {  
      for (const auto& spec : *pm_)
      {
        std::visit(overload {[&](const RangeStatsInt& /*int_range*/) { updateHistFromDataArray(spec.getIntegerDataArrays(), which.name, result); },
                             [&](const RangeStatsDouble& /*double_range*/) { updateHistFromDataArray(spec.getFloatDataArrays(), which.name, result); }}
                   , overview_range_data_.at(which));
      }
    }
    return result;
  }

  LayerStatisticsFeatureMap::LayerStatisticsFeatureMap(const FeatureMap& fm) : fm_(&fm)
  {
    computeStatistics_();
  }

  void addMetaDistributionValue(Math::Histogram<>& result, const string& name, const MetaInfoInterface& mi)
  {
    if (mi.metaValueExists(name))
    {
      result.inc(mi.getMetaValue(name));
    }
  }

  Math::Histogram<> LayerStatisticsFeatureMap::getDistribution(const RangeStatsType& which,
                                                               const UInt number_of_bins) const
  {
    auto mm = getMinMax(overview_range_data_, which, "FeatureMap"); // may throw if unknown statistic
    Math::Histogram<> result(mm.min, mm.max, (mm.max-mm.min) / number_of_bins);

    if (which.src == RangeStatsSource::CORE)
    {
      if (which.name == "intensity") for (const auto& f : *fm_) result.inc(f.getIntensity());
      else if (which.name == "charge") for (const auto& f : *fm_) result.inc(f.getCharge());
      else if (which.name == "quality") for (const auto& f : *fm_) result.inc(f.getOverallQuality());
    }
    else if (which.src == RangeStatsSource::METAINFO)
    {
      for (const auto& f : *fm_) addMetaDistributionValue(result, which.name, f);
    }

    return result;
  }
  void LayerStatisticsFeatureMap::computeStatistics_()
  {
    RangeStats<double> stat_intensity;
    RangeStats<int> stat_charge;
    RangeStats<double> stat_quality;
    for (const auto& f : *fm_)
    {
      stat_intensity.addDataPoint(f.getIntensity());
      stat_charge.addDataPoint(f.getCharge());
      stat_quality.addDataPoint(f.getOverallQuality());
      bringInMetaStats_(&f);
    }
    overview_range_data_.emplace(RangeStatsType {RangeStatsSource::CORE, "intensity"}, stat_intensity);
    overview_range_data_.emplace(RangeStatsType {RangeStatsSource::CORE, "charge"}, stat_charge);
    overview_range_data_.emplace(RangeStatsType {RangeStatsSource::CORE, "quality"}, stat_quality);
  }


  LayerStatisticsConsensusMap::LayerStatisticsConsensusMap(const ConsensusMap& cm) : cm_(&cm)
  {
    computeStatistics_();
  }

  Math::Histogram<> LayerStatisticsConsensusMap::getDistribution(const RangeStatsType& which,
                                                                 const UInt number_of_bins) const
  {
    auto mm = getMinMax(overview_range_data_, which, "ConsensusMap"); // may throw if unknown statistic
    Math::Histogram<> result(mm.min, mm.max, (mm.max - mm.min) / number_of_bins);

    if (which.src == RangeStatsSource::CORE)
    {
      if (which.name == "intensity") for (const auto& cf : *cm_) result.inc(cf.getIntensity());
      else if (which.name == "charge") for (const auto& cf : *cm_) result.inc(cf.getCharge());
      else if (which.name == "quality") for (const auto& cf : *cm_) result.inc(cf.getQuality());
      else if (which.name == "sub-elements") for (const auto& cf : *cm_) result.inc(cf.size());
    }
    else if (which.src == RangeStatsSource::METAINFO)
    {
      for (const auto& f : *cm_) addMetaDistributionValue(result, which.name, f);
    }
    
    return result;
  }

  void LayerStatisticsConsensusMap::computeStatistics_()
  {
    RangeStats<double> stat_intensity;
    RangeStats<int> stat_charge;
    RangeStats<double> stat_quality;
    RangeStats<int> stat_elements;
    for (const auto& cf : *cm_)
    {
      stat_intensity.addDataPoint(cf.getIntensity());
      stat_charge.addDataPoint(cf.getCharge());
      stat_quality.addDataPoint(cf.getQuality());
      stat_elements.addDataPoint(cf.size());
      bringInMetaStats_(&cf);
    }
    overview_range_data_.emplace(RangeStatsType {RangeStatsSource::CORE, "intensity"}, stat_intensity);
    overview_range_data_.emplace(RangeStatsType {RangeStatsSource::CORE, "charge"}, stat_charge);
    overview_range_data_.emplace(RangeStatsType {RangeStatsSource::CORE, "quality"}, stat_quality);
    overview_range_data_.emplace(RangeStatsType {RangeStatsSource::CORE, "sub-elements"},
                                 stat_quality);
  }

  // IDENT

  LayerStatisticsIdent::LayerStatisticsIdent(const IPeptideIds::PepIds& ids) : ids_(&ids)
  {
    computeStatistics_();
  }

  Math::Histogram<> LayerStatisticsIdent::getDistribution(const RangeStatsType& which,
                                                          const UInt number_of_bins) const
  {
    auto mm =
      getMinMax(overview_range_data_, which, "vector<PepIDs>"); // may throw if unknown statistic
    Math::Histogram<> result(mm.min, mm.max, (mm.max - mm.min) / number_of_bins);

    if (which.src == RangeStatsSource::METAINFO)
    {
      for (const auto& pep : *ids_)
        addMetaDistributionValue(result, which.name, pep);
    }

    return result;
  }

  void LayerStatisticsIdent::computeStatistics_()
  {
    for (const auto& pep : *ids_)
    {
      bringInMetaStats_(&pep);
    }
  }

  void LayerStatistics::bringInMetaStats_(const MetaInfoInterface* meta_interface)
  {
    vector<String> new_meta_keys;
    meta_interface->getKeys(new_meta_keys);
    for (const auto& idx : new_meta_keys)
    {
      const DataValue& meta_dv = meta_interface->getMetaValue(idx);
      if (meta_dv.valueType() == DataValue::INT_VALUE ||
          meta_dv.valueType() == DataValue::DOUBLE_VALUE)
      { // range data
        RangeStatsType key = {RangeStatsSource::METAINFO, idx};
        auto itr = overview_range_data_.find(key);
        // create correct type, if not present
        if (itr == overview_range_data_.end())
        {
          auto empty_value = [&]() -> RangeStatsVariant {
            if (meta_dv.valueType() == DataValue::INT_VALUE)
            {
              return RangeStatsInt();
            }
            else if (meta_dv.valueType() == DataValue::DOUBLE_VALUE)
            {
              return RangeStatsDouble();
            }
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Metavalue has unsupported valuetype", "??");}();// immediately evaluated lambda
          itr = overview_range_data_.emplace(key, empty_value).first; 
        }
        // update the value
        std::visit(overload {[&](auto&& stats) { stats.addDataPoint(meta_dv); }}, itr->second);        
      }
      else
      { // just count data
        ++overview_count_data_[idx].counter;
      }
    }
  }


} // namespace OpenMS