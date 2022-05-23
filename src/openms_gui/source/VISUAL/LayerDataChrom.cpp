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

#include <OpenMS/VISUAL/LayerDataChrom.h>

#include <OpenMS/VISUAL/LayerData1DPeak.h>
#include <OpenMS/VISUAL/LayerData1DChrom.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

#include <OpenMS/KERNEL/DimMapper.h>

using namespace std;

namespace OpenMS
{
  std::unique_ptr<LayerData1DBase> LayerDataChrom::to1DLayer() const
  {
    return make_unique<LayerData1DChrom>(*this);
  }

  std::unique_ptr<LayerStoreData> LayerDataChrom::storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const
  {
    auto ret = std::unique_ptr<LayerStoreDataPeakMapVisible>();
    ret->storeVisibleExperiment(*chromatogram_map_.get(), visible_range, layer_filters);
    return ret;
  }

  std::unique_ptr<LayerStoreData> LayerDataChrom::storeFullData() const
  {
    auto ret = std::unique_ptr<LayerStoreDataPeakMapAll>();
    ret->storeFullExperiment(*chromatogram_map_.get());
    return ret;
  }

  LayerDataChrom::ProjectionData LayerDataChrom::getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& area) const
  {
    ProjectionData result;

    // create projection data
    map<float, float> rt;
    map<int, float> mzint;
    map<int, float> mzsum;

    auto& peak_count = result.stats.number_of_datapoints;
    auto& intensity_max = result.stats.max_intensity;
    double total_intensity_sum = 0.0;

    // divide visible range into 100 bins (much faster than using a constant, e.g. 0.05, leading to many peaks for large maps without more information)
    float rt_range = area.RangeRT::getSpan();
    float mult = 100.0f / (std::isnan(rt_range) ? 1 : rt_range);

    MSSpectrum projection_mz;
    MSChromatogram projection_rt;

    // this does not work yet for chromatograms...
    /*for (auto i = chromatogram_map_->areaBeginConst(area.getMinRT(), area.getMaxRT(), area.getMinMZ(), area.getMaxMZ()); i != getPeakData()->areaEndConst(); ++i)
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
    ptr_rt->setPeakData(ExperimentSharedPtrType(new ExperimentType(exp_rt)));

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

  LayerDataDefs::PointXYType LayerDataChrom::peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const
  {
    return mapper.map(getChromatogram(peak.spectrum)[peak.peak]);
  }

  std::unique_ptr<LayerStatistics> LayerDataChrom::getStats() const
  {
    return make_unique<LayerStatisticsPeakMap>(*chromatogram_map_);
  }
} // namespace OpenMS
