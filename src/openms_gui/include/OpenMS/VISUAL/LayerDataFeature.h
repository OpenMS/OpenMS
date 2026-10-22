// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerDataBase.h>

namespace OpenMS
{

  /**
  @brief Class that stores the data for one layer of type FeatureMap

  The identifications of the features (and the unassigned ones) are the identification data of the map.

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataFeature : public virtual LayerDataBase
  {
  public:
    /// Default constructor
    LayerDataFeature();
    /// no Copy-ctor (should not be needed)
    LayerDataFeature(const LayerDataFeature& ld) = delete;
    /// no assignment operator (should not be needed)
    LayerDataFeature& operator=(const LayerDataFeature& ld) = delete;

    std::unique_ptr<Painter2DBase> getPainter2D() const override;

    std::unique_ptr<LayerData1DBase> to1DLayer() const override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }


    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;

    std::unique_ptr<LayerStoreData> storeFullData() const override;

    ProjectionData getProjection(const DIM_UNIT /*unit_x*/, const DIM_UNIT /*unit_y*/, const RangeAllType& /*area*/) const override
    { // currently only a stub
      ProjectionData proj;
      return proj;
    }

    PeakIndex findHighestDataPoint(const RangeAllType& area) const override;

    PointXYType peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const override;

    void updateRanges() override
    {
      features_->updateRanges();
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      r.assign(*getFeatureMap());
      return r;
    }

    std::unique_ptr<LayerStatistics> getStats() const override;

    bool annotate(const IdentificationData& identifications) override;


    /// Returns a const reference to the current feature data
    const FeatureMapSharedPtrType& getFeatureMap() const
    {
      return features_;
    }

    /// Returns a const reference to the current feature data
    FeatureMapSharedPtrType& getFeatureMap()
    {
      return features_;
    }

  protected:
    /// feature data
    FeatureMapSharedPtrType features_ = FeatureMapSharedPtrType(new FeatureMapType());
  };

}// namespace OpenMS
