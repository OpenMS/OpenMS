// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerDataBase.h>

#include <vector>

namespace OpenMS
{

  /**
  @brief Class that stores the data for one layer of type ConsensusMap

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataConsensus : public virtual LayerDataBase
  {
  public:
    /// Default constructor
    LayerDataConsensus(ConsensusMapSharedPtrType& map);
    /// no Copy-ctor (should not be needed)
    LayerDataConsensus(const LayerDataConsensus& ld) = delete;
    /// no assignment operator (should not be needed)
    LayerDataConsensus& operator=(const LayerDataConsensus& ld) = delete;

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
      consensus_map_->updateRanges();
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      r.assign(*getConsensusMap());
      return r;
    }

    std::unique_ptr<LayerStatistics> getStats() const override;

    bool annotate(const std::vector<PeptideIdentification>& identifications, const std::vector<ProteinIdentification>& protein_identifications) override;
  
    /// Returns a const reference to the consensus feature data
    const ConsensusMapSharedPtrType& getConsensusMap() const
    {
      return consensus_map_;
    }

    /// Returns current consensus map (mutable)
    ConsensusMapSharedPtrType& getConsensusMap()
    {
      return consensus_map_;
    }
  protected:
    /// consensus feature data
    ConsensusMapSharedPtrType consensus_map_ = ConsensusMapSharedPtrType(new ConsensusMapType());
  };

}// namespace OpenMS
