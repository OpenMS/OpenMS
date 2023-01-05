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

#include <OpenMS/KERNEL/Mobilogram.h>

namespace OpenMS
{
  
  /**
  @brief Class that stores the data for one layer of type IonMobility

  FIXME: currently we only store a single mobilogram, since this is what is required to show a projection in 2D View.
  If there is another application, feel free to implement a surrounding container (don't use vector<Mobilogram>, do it properly! :)).

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataIonMobility : public virtual LayerDataBase
  {
  public:
    using PeakType = Mobilogram::PeakType;

    /// Default constructor
    LayerDataIonMobility();
    /// Copy-ctor
    LayerDataIonMobility(const LayerDataIonMobility& ld);
    /// Assignment operator
    LayerDataIonMobility& operator=(const LayerDataIonMobility& ld) = default;
    /// move Ctor
    LayerDataIonMobility(LayerDataIonMobility&& ld) = default;
    /// move assignment
    LayerDataIonMobility& operator=(LayerDataIonMobility&& ld) = default;

    std::unique_ptr<Painter2DBase> getPainter2D() const override;

    std::unique_ptr<LayerData1DBase> to1DLayer() const override;

    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;

    std::unique_ptr<LayerStoreData> storeFullData() const override;

    ProjectionData getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& area) const override;

    PeakIndex findHighestDataPoint(const RangeAllType& /*area*/) const override
    { // todo: not implemented
      return PeakIndex();
    }

    void updateRanges() override
    {
      single_mobilogram_.updateRanges();
      // on_disc_peaks->updateRanges(); // note: this is not going to work since its on disk! We currently don't have a good way to access these ranges
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      r.assign(single_mobilogram_);
      return r;
    }

    // for now, only a single Mobilogram. See class description.
    void setMobilityData(const Mobilogram& mobilogram)
    {
      single_mobilogram_ = mobilogram;
    }

    const Mobilogram& getMobilogram(Size index) const
    {
      if (index != 0) throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Only one mobilogram possible atm.", String(index));
      return single_mobilogram_;
    }

    PointXYType peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const override;

    String getDataArrayDescription(const PeakIndex& peak_index) override;

    std::unique_ptr<LayerStatistics> getStats() const override;

  protected:
    Mobilogram single_mobilogram_; ///< a single mobilogram (for now) -- see class description
  };

}// namespace OpenMS
