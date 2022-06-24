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

#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/LayerData1DBase.h>

namespace OpenMS
{
  class Annotation1DItem;

  /**
  @brief Class that stores the data for one layer of type PeakMap

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataPeak : public virtual LayerDataBase
  {
  public:

    using SpectrumType = ExperimentType::SpectrumType;
    using PeakType = SpectrumType::PeakType;

    /// Default constructor
    LayerDataPeak();
    /// Copy-ctor
    LayerDataPeak(const LayerDataPeak& ld);
    /// Assignment operator
    LayerDataPeak& operator=(const LayerDataPeak& ld) = default;
    /// move Ctor
    LayerDataPeak(LayerDataPeak&& ld) = default;
    /// move assignment
    LayerDataPeak& operator=(LayerDataPeak&& ld) = default;

    std::unique_ptr<LayerData1DBase> to1DLayer() const override;

    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;

    std::unique_ptr<LayerStoreData> storeFullData() const override;

    ProjectionData getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& area) const override;

    void updateRanges() override
    {
      peak_map_->updateRanges();
      // on_disc_peaks->updateRanges(); // note: this is not going to work since its on disk! We currently don't have a good way to access these ranges
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      r.assign(*peak_map_);
      return r;
    }

    const ExperimentType::SpectrumType& getSpectrum(Size spectrum_idx) const
    {
      if ((*peak_map_)[spectrum_idx].size() > 0)
      {
        return (*peak_map_)[spectrum_idx];
      }
      if (!on_disc_peaks->empty())
      {
        static MSSpectrum local_spec;
        local_spec = on_disc_peaks->getSpectrum(spectrum_idx);
        return local_spec;
      }
      return (*peak_map_)[spectrum_idx];
    }

    PointXYType peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const override;

    String getDataArrayDescription(const PeakIndex& peak_index) override;

    std::unique_ptr<LayerStatistics> getStats() const override;
  };

}// namespace OpenMS
