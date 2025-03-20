// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerData1DBase.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>

namespace OpenMS
{
  
  class OPENMS_GUI_DLLAPI LayerData1DPeak : public LayerDataPeak, public LayerData1DBase
  {
  public:
    LayerData1DPeak()
      : LayerDataBase(DT_PEAK)
    {
    }

    LayerData1DPeak(const LayerDataPeak& base) : LayerDataBase(base), LayerDataPeak(base)
    {
    }


    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;
    std::unique_ptr<LayerStoreData> storeFullData() const override;

    std::unique_ptr<Painter1DBase> getPainter1D() const override;

    bool hasIndex(Size index) const override
    {
      return index < peak_map_->size();
    }

    RangeAllType getRangeForArea(const RangeAllType partial_range) const override
    {
      const auto& spec = getCurrentSpectrum();
      auto spec_filtered = SpectrumType();
      spec_filtered.insert(spec_filtered.begin(), spec.MZBegin(partial_range.getMinMZ()), spec.MZEnd(partial_range.getMaxMZ()));
      spec_filtered.updateRanges();
      return RangeAllType().assign(spec_filtered.getRange());
    }

    const ExperimentType::SpectrumType& getCurrentSpectrum() const
    {
      return LayerDataPeak::getSpectrum(current_idx_);
    }


    void updateRanges() override
    {
      LayerDataPeak::updateRanges();
    }

    RangeAllType getRange1D() const override
    {
      return RangeAllType().assign(getCurrentSpectrum().getRange());
    }

    // docu in base class
    QMenu* getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint) override;

    PeakIndex findClosestDataPoint(const RangeAllType& area) const override;

    // docu in base class
    Annotation1DItem* addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color) override;

    /// updates the PeakAnnotations in the current PeptideHit with manually changed annotations
    /// if no PeptideIdentification or PeptideHit for the spectrum exist, it is generated
    void synchronizePeakAnnotations();

    /// remove peak annotations in the given list from the currently active PeptideHit
    void removePeakAnnotationsFromPeptideHit(const std::vector<Annotation1DItem*>& selected_annotations);

    /// updates the PeakAnnotations in the current PeptideHit with manually changed annotations
    void updatePeptideHitAnnotations_(PeptideHit& hit);

  protected:
  };

}// namespace OpenMS
