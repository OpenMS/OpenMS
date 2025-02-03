// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/VISUAL/LayerData1DBase.h>
#include <OpenMS/VISUAL/LayerDataChrom.h>

namespace OpenMS
{
  class OPENMS_GUI_DLLAPI LayerData1DChrom : public LayerDataChrom, public LayerData1DBase
  {
  public:
    LayerData1DChrom() : LayerDataBase(DT_CHROMATOGRAM)
    {
    }

    LayerData1DChrom(const LayerDataChrom& base) : LayerDataBase(base), LayerDataChrom(base)
    {
    }

    std::unique_ptr<LayerStoreData> storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const override;
    std::unique_ptr<LayerStoreData> storeFullData() const override;

    std::unique_ptr<Painter1DBase> getPainter1D() const override;

    bool hasIndex(Size index) const override
    {
      return index < chromatogram_map_->getNrChromatograms();
    }

    RangeAllType getRangeForArea(const RangeAllType partial_range) const override
    {
      // update ranges based on given RT range
      if (partial_range.RangeRT::isEmpty())
      { // .. unless RT is empty, then we use the whole RT range
        auto r = RangeAllType(partial_range);
        r.extend(getCurrentChrom().getRange());
        return r;
      }
      const auto& chrom = getCurrentChrom();
      auto chrom_filtered = MSExperiment::ChromatogramType();
      chrom_filtered.insert(chrom_filtered.begin(), chrom.RTBegin(partial_range.getMinRT()), chrom.RTEnd(partial_range.getMaxRT()));
      chrom_filtered.updateRanges();
      return RangeAllType().assign(chrom_filtered.getRange());
    }

    RangeAllType getRange1D() const override
    {
      return RangeAllType().assign(getCurrentChrom().getRange());
    }

    const ExperimentType::ChromatogramType& getCurrentChrom() const
    {
      return getChromatogram(current_idx_);
    }

    void updateRanges() override
    {
      LayerDataChrom::updateRanges();
    }

    RangeAllType getRange() const override
    {
      // do NOT change the behaviour of getRange() for 1D, since we want the full RT range across all chroms
      // when scrolling in the list of chroms
      return LayerDataChrom::getRange();
    }

    // docu in base class
    QMenu* getContextMenuAnnotation(Annotation1DItem* annot_item, bool& need_repaint) override;

    PeakIndex findClosestDataPoint(const RangeAllType& area) const override;

    // docu in base class
    Annotation1DItem* addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color) override;

  protected:
    /// Current cached spectrum
    //ExperimentType::SpectrumType cached_spectrum_;
  };

} // namespace OpenMS
