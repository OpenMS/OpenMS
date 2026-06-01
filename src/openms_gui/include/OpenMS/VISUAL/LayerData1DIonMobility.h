// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerData1DBase.h>
#include <OpenMS/VISUAL/LayerDataIonMobility.h>

namespace OpenMS
{
  
  class OPENMS_GUI_DLLAPI LayerData1DIonMobility : public LayerDataIonMobility, public LayerData1DBase
  {
  public:
    LayerData1DIonMobility()
      : LayerDataBase(DT_PEAK)
    {
    }

    LayerData1DIonMobility(const LayerDataIonMobility& base) : LayerDataBase(base), LayerDataIonMobility(base)
    {
    }


    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;
    std::unique_ptr<LayerStoreData> storeFullData() const override;

    std::unique_ptr<Painter1DBase> getPainter1D() const override;

    bool hasIndex(Size index) const override
    {
      return index == 0;
    }

    RangeAllType getRangeForArea(const RangeAllType partial_range) const override
    {
      const auto& spec = getCurrentMobilogram();
      auto spec_filtered = Mobilogram();
      spec_filtered.insert(spec_filtered.begin(), spec.MBBegin(partial_range.getMinMobility()), spec.MBEnd(partial_range.getMaxMobility()));
      spec_filtered.updateRanges();
      return RangeAllType().assign(spec_filtered.getRange());
    }
    
    RangeAllType getRange1D() const override
    {
      return RangeAllType().assign(getCurrentMobilogram().getRange());
    }

    const Mobilogram& getCurrentMobilogram() const
    {
      return LayerDataIonMobility::getMobilogram(this->getCurrentIndex());
    }

    void updateRanges() override
    {
      LayerDataIonMobility::updateRanges();
    }

    RangeAllType getRange() const override
    {
      // do NOT change the behaviour of getRange() for 1D, since we want the full IM range across all mbs
      // when scrolling in the list of mbs
      return LayerDataIonMobility::getRange();
    }

    // docu in base class
    QMenu* getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint) override;

    PeakIndex findClosestDataPoint(const RangeAllType& area) const override;

    // docu in base class
    Annotation1DItem* addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color) override;

  protected:
  };

}// namespace OpenMS
