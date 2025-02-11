// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/INTERFACES/IPeptideIds.h>

namespace OpenMS
{

  /**
  @brief Class that stores the data for one layer of type PeptideIdentifications

  @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI LayerDataIdent : public LayerDataBase, public IPeptideIds
  {
  public:
    /// Default constructor
    LayerDataIdent() :
        LayerDataBase(LayerDataBase::DT_IDENT){};
    /// no Copy-ctor (should not be needed)
    LayerDataIdent(const LayerDataIdent& ld) = delete;
    /// no assignment operator (should not be needed)
    LayerDataIdent& operator=(const LayerDataIdent& ld) = delete;

    std::unique_ptr<Painter2DBase> getPainter2D() const override;

    std::unique_ptr<LayerData1DBase> to1DLayer() const override
    {
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;

    std::unique_ptr<LayerStoreData> storeFullData() const override;

    ProjectionData getProjection(const DIM_UNIT unit_x, const DIM_UNIT unit_y, const RangeAllType& area) const override;

    PeakIndex findHighestDataPoint(const RangeAllType& /*area*/) const override
    { // todo: not implemented
      return PeakIndex();
    }

    PointXYType peakIndexToXY(const PeakIndex& peak, const DimMapper<2>& mapper) const override;


    void updateRanges() override
    {
      // nothing to do...
    }

    RangeAllType getRange() const override
    {
      RangeAllType r;
      for (const PeptideIdentification& pep : peptides_)
      {
        r.extendRT(pep.getRT());
        r.extendMZ(pep.getMZ());
      }
      return r;
    }

    std::unique_ptr<LayerStatistics> getStats() const override;

    virtual const PepIds& getPeptideIds() const override
    {
      return peptides_;
    }
    virtual PepIds& getPeptideIds() override
    {
      return peptides_;
    }

    virtual void setPeptideIds(const PepIds& ids) override
    {
      peptides_ = ids;
    }
    virtual void setPeptideIds(PepIds&& ids) override
    {
      peptides_ = std::move(ids);
    }

  private:
    /// peptide identifications
    std::vector<PeptideIdentification> peptides_;
  };

}// namespace OpenMS
