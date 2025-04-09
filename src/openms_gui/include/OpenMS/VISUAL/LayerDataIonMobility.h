// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
    /// no assignment operator (should not be needed)
    LayerDataIonMobility& operator=(const LayerDataIonMobility& ld) = delete;


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
