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

// OpenMS includes
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/MetaInfo.h>

using namespace std;

namespace OpenMS
{

  /// Computes the statistics of all meta data contained in the FloatDataArray or IntegerDataArray of an MSSpectrum
  template<typename SSType, typename DataArrayType>
  void computeMetaDataArrayStats_(const DataArrayType& arrays, std::map<std::string, StatsSummaryVariant>& meta_array_stats_)
  {
    using VarType = StatsSummary<SSType>;

    for (const auto& mda : arrays)
    {
      const String& mda_name = mda.getName();
      auto it_var = meta_array_stats_.find(mda_name);
      if (it_var == meta_array_stats_.end()) // missing, create it
      {
        it_var = meta_array_stats_.emplace(mda_name, VarType()).first;
      }
      StatsSummaryVariant& meta_stats_value = it_var->second;
      auto& target = std::get<VarType>(meta_stats_value);
      for (const auto& value : mda)
      {
        target.addDataPoint(value);
      }
    }
  }

  void LayerStatistics::computePeakMapStats(const PeakMap& pm)
  {
    StatsSummary<double> stat_intensity;
    for (const auto& spec : pm)
    {
      for (const auto& peak : spec)
      {
        stat_intensity.addDataPoint(peak.getIntensity());
      }
      // collect stats about the meta data arrays of this spectrum
      computeMetaDataArrayStats_<double>(spec.getFloatDataArrays(), meta_array_stats_);
      computeMetaDataArrayStats_<int>(spec.getIntegerDataArrays(), meta_array_stats_);
    }
    core_stats_.emplace("intensity", stat_intensity);
  }

  void LayerStatistics::computeFeatureMapStats(const FeatureMap& fm)
  {
    StatsSummary<double> stat_intensity;
    StatsSummary<int> stat_charge;
    StatsSummary<double> stat_quality;
    for (const auto& f : fm)
    {
      stat_charge.addDataPoint(f.getCharge());
      stat_quality.addDataPoint(f.getOverallQuality());
      stat_intensity.addDataPoint(f.getIntensity());
      bringInMetaStats_(&f);
    }
    core_stats_.emplace("intensity", stat_intensity);
    core_stats_.emplace("charge", stat_charge);
    core_stats_.emplace("quality", stat_quality);
  }

  void LayerStatistics::computeConsensusMapStats(const ConsensusMap& cm)
  {
    StatsSummary<double> stat_intensity;
    StatsSummary<int> stat_charge;
    StatsSummary<int> stat_elements;
    StatsSummary<double> stat_quality;
    for (const auto& cf : cm)
    {
      stat_charge.addDataPoint(cf.getCharge());
      stat_quality.addDataPoint(cf.getQuality());
      stat_intensity.addDataPoint(cf.getIntensity());
      stat_elements.addDataPoint(cf.size());
      bringInMetaStats_(&cf);
    }
    core_stats_.emplace("intensity", stat_intensity);
    core_stats_.emplace("charge", stat_charge);
    core_stats_.emplace("quality", stat_quality);
    core_stats_.emplace("sub-elements", stat_quality);
  }

  /// retrieve the core statistics of the data structure (usually intensity, maybe charge, quality, etc)

  const LayerStatistics::StatsMap& LayerStatistics::getCoreStats() const
  {
    return core_stats_;
  }

  /// retrieve the statistics of all the meta values, e.g. "FWHM" for FeatureMaps
  /// The result is computed by looking up a MetaInfoInterface. Try not to query repeatedly if performance matters.

  LayerStatistics::StatsMap LayerStatistics::getMetaStats() const
  {
    // translate from UInt to string
    StatsMap result;
    for (const auto& entry : meta_stats_)
    {
      result.emplace(MetaInfo::registry().getName(entry.first), entry.second);
    }
    return result;
  }

  /// retrieve the statistics of all float and integer data arrays of spectra (only populated for PeakMaps)

  const LayerStatistics::StatsMap& LayerStatistics::getArrayStats() const
  {
    return meta_array_stats_;
  }

  // helper for visitor pattern with std::visit
  template<class... Ts>
  struct overload : Ts... {
    using Ts::operator()...;
  };
  template<class... Ts>
  overload(Ts...) -> overload<Ts...>;


  void LayerStatistics::bringInMetaStats_(const MetaInfoInterface* meta_interface)
  {
    vector<UInt> new_meta_keys;
    meta_interface->getKeys(new_meta_keys);
    for (const UInt idx : new_meta_keys)
    {
      const DataValue& meta_dv = meta_interface->getMetaValue(idx);
      auto it = meta_stats_.find(idx);
      // create correct type, if not present
      if (it == meta_stats_.end())
      {
        StatsSummaryVariant v = [&]() -> StatsSummaryVariant
        {
          if (meta_dv.valueType() == DataValue::INT_VALUE)
          {
            return SSInt();
          }
          else if (meta_dv.valueType() == DataValue::DOUBLE_VALUE)
          {
            return SSDouble();
          }
          else
          { // simply count how often a metavalue occurs (e.g. a DataValue::String)
            return StatsCounter();
          }
        }();
        it = meta_stats_.emplace(idx, v).first;
      }
      auto& val = it->second;
      // update the value
      std::visit(overload{
            [&](auto&& stats) { stats.addDataPoint(meta_dv); },
            [&](StatsCounter& c) { ++c.counter; }
      }, val);
    }
  }

} // namespace OpenMS