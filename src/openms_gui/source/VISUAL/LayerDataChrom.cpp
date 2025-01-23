// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------



#include <OpenMS/VISUAL/LayerDataChrom.h>

#include <OpenMS/DATASTRUCTURES/OSWData.h>
#include <OpenMS/KERNEL/DimMapper.h>

#include <OpenMS/VISUAL/LayerData1DPeak.h>
#include <OpenMS/VISUAL/LayerData1DChrom.h>
#include <OpenMS/VISUAL/Painter2DBase.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>


using namespace std;

namespace OpenMS
{
  LayerDataChrom::LayerDataChrom():
    LayerDataBase(LayerDataBase::DT_CHROMATOGRAM)
  {}

  std::unique_ptr<Painter2DBase> LayerDataChrom::getPainter2D() const
  {
    return make_unique<Painter2DChrom>(this);
  }

  std::unique_ptr<LayerData1DBase> LayerDataChrom::to1DLayer() const
  {
    return make_unique<LayerData1DChrom>(*this);
  }

  std::unique_ptr<LayerStoreData> LayerDataChrom::storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const
  {
    auto ret = make_unique<LayerStoreDataPeakMapVisible>();
    ret->storeVisibleExperiment(*chromatogram_map_.get(), visible_range, layer_filters);
    return ret;
  }

  std::unique_ptr<LayerStoreData> LayerDataChrom::storeFullData() const
  {
    auto ret = make_unique<LayerStoreDataPeakMapAll>();
    ret->storeFullExperiment(*chromatogram_map_.get());
    return ret;
  }

  LayerDataChrom::ProjectionData LayerDataChrom::getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& /*area*/) const
  {
    ProjectionData result;

    // create projection data
    map<float, float> rt;
    map<int, float> mzint;
    map<int, float> mzsum;

    MSSpectrum projection_mz;
    MSChromatogram projection_rt;
        
    // this does not work yet for chromatograms...
    /*
    auto& peak_count = result.stats.number_of_datapoints;
    auto& intensity_max = result.stats.max_intensity;
    double total_intensity_sum = 0.0;

    // divide visible range into 100 bins (much faster than using a constant, e.g. 0.05, leading to many peaks for large maps without more information)
    float rt_range = area.RangeRT::getSpan();
    float mult = 100.0f / (std::isnan(rt_range) ? 1 : rt_range);

    for (auto i = chromatogram_map_->areaBeginConst(area.getMinRT(), area.getMaxRT(), area.getMinMZ(), area.getMaxMZ()); i != getPeakData()->areaEndConst(); ++i)
    {
      PeakIndex pi = i.getPeakIndex();
      if (filters.passes((*getPeakData())[pi.spectrum], pi.peak))
      {
        // summary stats
        ++peak_count;
        total_intensity_sum += i->getIntensity();
        intensity_max = max(intensity_max, i->getIntensity());

        // binning for m/z
        auto intensity = i->getIntensity();
        mzint[int(i->getMZ() * mult)] += intensity;
        // ... to later obtain an intensity weighted average m/z value
        mzsum[int(i->getMZ() * mult)] += i->getMZ() * intensity;

        // binning in RT (one value per scan)
        rt[i.getRT()] += i->getIntensity();
      }
    }

    // write to spectra/chrom
    projection_mz.resize(mzint.size() + 2);
    projection_mz[0].setMZ(area.getMinMZ());
    projection_mz[0].setIntensity(0.0);
    projection_mz.back().setMZ(area.getMaxMZ());
    projection_mz.back().setIntensity(0.0);

    projection_rt.resize(rt.size() + 2);
    projection_rt[0].setRT(area.getMinRT());
    projection_rt[0].setIntensity(0.0);
    projection_rt.back().setRT(area.getMaxRT());
    projection_rt.back().setIntensity(0.0);

    Size i = 1;
    map<int, float>::iterator intit = mzint.begin();

    for (auto it = mzsum.cbegin(); it != mzsum.cend(); ++it)
    {
      auto intensity = intit->second;
      projection_mz[i].setMZ(it->second / intensity);
      projection_mz[i].setIntensity(intensity);
      ++intit;
      ++i;
    }

    i = 1;
    for (auto it = rt.cbegin(); it != rt.cend(); ++it)
    {
      projection_rt[i].setMZ(it->first);
      projection_rt[i].setIntensity(it->second);
      ++i;
    }
    */

    // create final datastructure

    // projection for m/z
    auto ptr_mz = make_unique<LayerData1DPeak>();
    MSExperiment exp_mz;
    exp_mz.addSpectrum(std::move(projection_mz));
    ptr_mz->setPeakData(ExperimentSharedPtrType(new ExperimentType(exp_mz)));

    // projection for RT
    auto ptr_rt = make_unique<LayerData1DChrom>();
    MSExperiment exp_rt;
    exp_mz.addChromatogram(std::move(projection_rt));
    ptr_rt->setChromData(ExperimentSharedPtrType(new ExperimentType(exp_rt)));

    auto assign_axis = [&](auto unit, auto& layer) {
      switch (unit)
      {
        case DIM_UNIT::MZ:
          layer = std::move(ptr_mz);
          break;
        case DIM_UNIT::RT:
          layer = std::move(ptr_rt);
          break;
        default:
          // do nothing, leave projection empty
          break;
      }
    };
    assign_axis(unit_x, result.projection_ontoX);
    assign_axis(unit_y, result.projection_ontoY);

    return result;
  }

  PeakIndex LayerDataChrom::findHighestDataPoint(const RangeAllType& area) const
  {
    const PeakMap& exp = *getChromatogramData();
    int count {-1};
    for (const auto& chrom : exp.getChromatograms())
    {
      ++count;
      if (chrom.empty())
      {
        continue; // ensure that empty chromatograms are not examined (iter->front = segfault)
      }

      auto mz_origin = chrom.getPrecursor().getMZ();

      // check m/z first
      if (area.containsMZ(mz_origin))
      {
        // the center point's RT should be inside the RT range of the chromatogram
        if (RangeRT(chrom.front().getRT(), chrom.back().getRT()).containsRT(area.RangeRT::center()))
        {
          return PeakIndex(count, 0); // we only care about the chrom, not the peak inside
        }
      }
    }
    return PeakIndex();
  }

  PointXYType LayerDataChrom::peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const
  {
    const auto& chrom = getChromatogram(peak.spectrum);
    return mapper.map(chrom, peak.peak);
  }

  String LayerDataChrom::getDataArrayDescription(const PeakIndex& peak_index)
  {
    String status;
    const auto& s = getChromatogram(peak_index.spectrum);
    for (Size m = 0; m < s.getFloatDataArrays().size(); ++m)
    {
      if (peak_index.peak < s.getFloatDataArrays()[m].size())
      {
        status += s.getFloatDataArrays()[m].getName() + ": " + s.getFloatDataArrays()[m][peak_index.peak] + " ";
      }
    }
    for (Size m = 0; m < s.getIntegerDataArrays().size(); ++m)
    {
      if (peak_index.peak < s.getIntegerDataArrays()[m].size())
      {
        status += s.getIntegerDataArrays()[m].getName() + ": " + s.getIntegerDataArrays()[m][peak_index.peak] + " ";
      }
    }
    for (Size m = 0; m < s.getStringDataArrays().size(); ++m)
    {
      if (peak_index.peak < s.getStringDataArrays()[m].size())
      {
        status += s.getStringDataArrays()[m].getName() + ": " + s.getStringDataArrays()[m][peak_index.peak] + " ";
      }
    }
    return status;
  }

  void LayerDataChrom::setChromatogramAnnotation(OSWData&& data)
  {
    chrom_annotation_ = OSWDataSharedPtrType(new OSWData(std::move(data)));
  }

  std::unique_ptr<LayerStatistics> LayerDataChrom::getStats() const
  {
    return make_unique<LayerStatisticsPeakMap>(*chromatogram_map_);
  }
} // namespace OpenMS
