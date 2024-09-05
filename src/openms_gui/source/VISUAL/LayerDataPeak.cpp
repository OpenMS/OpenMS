// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerDataPeak.h>

#include <OpenMS/ANALYSIS/ID/IDMapper.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/LayerData1DIonMobility.h>
#include <OpenMS/VISUAL/LayerData1DChrom.h>
#include <OpenMS/VISUAL/LayerData1DPeak.h>
#include <OpenMS/VISUAL/Painter2DBase.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

using namespace std;

namespace OpenMS
{
  LayerDataPeak::LayerDataPeak() : LayerDataBase(LayerDataBase::DT_PEAK)
  {
    flags.set(LayerDataBase::P_PRECURSORS);
  }

  /*LayerDataPeak::LayerDataPeak(const LayerDataPeak& ld)
    : LayerDataBase(static_cast<const LayerDataBase&>(ld)),
      peak_map_(ld.peak_map_),
      on_disc_peaks_(ld.on_disc_peaks_)
  {
  } */

  std::unique_ptr<Painter2DBase> LayerDataPeak::getPainter2D() const
  {
    return make_unique<Painter2DPeak>(this);
  }

  std::unique_ptr<LayerData1DBase> LayerDataPeak::to1DLayer() const
  {
    return make_unique<LayerData1DPeak>(*this);
  }

  std::unique_ptr<LayerStoreData> LayerDataPeak::storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const
  {
    auto ret = make_unique<LayerStoreDataPeakMapVisible>();
    ret->storeVisibleExperiment(*peak_map_.get(), visible_range, layer_filters);
    return ret;
  }

  std::unique_ptr<LayerStoreData> LayerDataPeak::storeFullData() const
  {
    auto ret = make_unique<LayerStoreDataPeakMapAll>();
    ret->storeFullExperiment(*peak_map_.get());
    return ret;
  }

  LayerDataPeak::ProjectionData LayerDataPeak::getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& area) const
  {
    ProjectionData result;

    // create projection data
    map<float, float> rt;
    map<float, float> mobility;
    map<int, float> mzint;
    map<int, float> mzsum;

    auto& peak_count = result.stats.number_of_datapoints;
    auto& intensity_max = result.stats.max_intensity;
    auto& total_intensity_sum = result.stats.sum_intensity; 

    // divide visible range into 100 bins (much faster than using a constant, e.g. 0.05, leading to many peaks for large maps without more information)
    float mz_range = area.RangeMZ::getSpan();
    float mult = 100.0f / (std::isnan(mz_range) ? 1 : mz_range);

    MSSpectrum projection_mz;
    Mobilogram projection_im;
    MSChromatogram projection_rt;

    for (auto i = getPeakData()->areaBeginConst(area); i != getPeakData()->areaEndConst(); ++i)
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

        // binning in IM (one value per scan)
        mobility[i.getDriftTime()] += i->getIntensity();
      }
    }

    // write to spectra/chrom
    projection_mz.resize(mzint.size() + 2);
    projection_mz[0].setMZ(area.getMinMZ());
    projection_mz[0].setIntensity(0.0);
    projection_mz.back().setMZ(area.getMaxMZ());
    projection_mz.back().setIntensity(0.0);

    projection_im.resize(mobility.size() + 2);
    projection_im[0].setMobility(area.getMinMobility());
    projection_im[0].setIntensity(0.0);
    projection_im.back().setMobility(area.getMaxMobility());
    projection_im.back().setIntensity(0.0);
    

    projection_rt.resize(rt.size() + 2);
    projection_rt[0].setRT(area.getMinRT());
    projection_rt[0].setIntensity(0.0);
    projection_rt.back().setRT(area.getMaxRT());
    projection_rt.back().setIntensity(0.0);

    Size i = 1;
    auto intit = mzint.begin();
    for (auto it = mzsum.cbegin(); it != mzsum.cend(); ++it)
    {
      auto intensity = intit->second;
      projection_mz[i].setMZ(it->second / intensity);
      projection_mz[i].setIntensity(intensity);
      ++intit;
      ++i;
    }

    i = 1;
    for (auto it = mobility.cbegin(); it != mobility.cend(); ++it)
    {
      projection_im[i].setMobility(it->first);
      projection_im[i].setIntensity(it->second);
      ++i;
    }


    i = 1;
    for (auto it = rt.cbegin(); it != rt.cend(); ++it)
    {
      projection_rt[i].setMZ(it->first);
      projection_rt[i].setIntensity(it->second);
      ++i;
    }

    // create final datastructure

    // projection for m/z
    auto ptr_mz = make_unique<LayerData1DPeak>();
    {
      MSExperiment exp_mz;
      exp_mz.addSpectrum(std::move(projection_mz));
      ptr_mz->setPeakData(ExperimentSharedPtrType(new ExperimentType(std::move(exp_mz))));
    }

    // projection for mobility
    auto ptr_im = make_unique<LayerData1DIonMobility>();
    {
      ptr_im->setMobilityData(projection_im);
    }

    // projection for RT
    auto ptr_rt = make_unique<LayerData1DChrom>();
    {
      MSExperiment exp_rt;
      exp_rt.addChromatogram(std::move(projection_rt));
      ptr_rt->setChromData(ExperimentSharedPtrType(new ExperimentType(std::move(exp_rt))));
    }

    auto assign_axis = [&](auto unit, auto& layer) {
      switch (unit)
      {
        case DIM_UNIT::MZ:
          layer = std::move(ptr_mz);
          break;

        case DIM_UNIT::RT:
          layer = std::move(ptr_rt);
          break;

        case DIM_UNIT::FAIMS_CV:
        case DIM_UNIT::IM_MS:
        case DIM_UNIT::IM_VSSC:
          layer = std::move(ptr_im);
          break;

        default:
          throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
    };
    assign_axis(unit_x, result.projection_ontoX);
    assign_axis(unit_y, result.projection_ontoY);

    return result;
  }

  PeakIndex LayerDataPeak::findHighestDataPoint(const RangeAllType& area) const
  {
    using IntType = MSExperiment::ConstAreaIterator::PeakType::IntensityType;
    auto max_int = numeric_limits<IntType>::lowest();
    PeakIndex max_pi;
    
    const auto map = *getPeakData();
    // for IM data, use whatever is there. For RT/mz data, use MSlevel 1
    const UInt MS_LEVEL = (! map.empty() && map.isIMFrame()) ? map[0].getMSLevel() : 1;

    for (ExperimentType::ConstAreaIterator i = map.areaBeginConst(area, MS_LEVEL); i != map.areaEndConst(); ++i)
    {
      PeakIndex pi = i.getPeakIndex();
      if (i->getIntensity() > max_int && filters.passes((map)[pi.spectrum], pi.peak))
      {
        max_int = i->getIntensity();
        max_pi = pi;
      }
    }
    return max_pi;
  }

  PointXYType LayerDataPeak::peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const
  {
    const auto& spec = getSpectrum(peak.spectrum);
    return mapper.map(spec, peak.peak);
  }

  String LayerDataPeak::getDataArrayDescription(const PeakIndex& peak_index)
  {
    String status;
    const ExperimentType::SpectrumType& s = getSpectrum(peak_index.spectrum);
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

  std::unique_ptr<LayerStatistics> LayerDataPeak::getStats() const
  {
    return make_unique<LayerStatisticsPeakMap>(*peak_map_);
  }

  bool LayerDataPeak::annotate(const vector<PeptideIdentification>& identifications, const vector<ProteinIdentification>& protein_identifications)
  {
    IDMapper mapper;
    Param p = mapper.getDefaults();
    p.setValue("rt_tolerance", 0.1, "RT tolerance (in seconds) for the matching");
    p.setValue("mz_tolerance", 1.0, "m/z tolerance (in ppm or Da) for the matching");
    p.setValue("mz_measure", "Da", "unit of 'mz_tolerance' (ppm or Da)");
    mapper.setParameters(p);
    mapper.annotate(*getPeakDataMuteable(), identifications, protein_identifications, true);

    return true;
  }
  
  const LayerDataBase::ConstExperimentSharedPtrType LayerDataPeak::getPeakData() const
  {
    return boost::static_pointer_cast<const ExperimentType>(peak_map_);
  }

} // namespace OpenMS