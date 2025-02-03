// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/LayerData1DChrom.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/Painter1DBase.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>
#include <QMenu>
using namespace std;

namespace OpenMS
{
  std::unique_ptr<LayerStoreData> LayerData1DChrom::storeVisibleData(const RangeAllType& visible_range, const DataFilters& layer_filters) const
  {
    auto ret = std::make_unique<LayerStoreDataPeakMapVisible>();
    ret->storeVisibleChromatogram(getCurrentChrom(), visible_range, layer_filters);
    return ret;
  }

  std::unique_ptr<LayerStoreData> LayerData1DChrom::storeFullData() const
  {
    return LayerDataChrom::storeFullData(); // just forward
  }

  QMenu* LayerData1DChrom::getContextMenuAnnotation(Annotation1DItem* /*annot_item*/, bool& /*need_repaint*/)
  {
    auto* context_menu = new QMenu("Chrom1D", nullptr);

    return context_menu;
  }

  PeakIndex LayerData1DChrom::findClosestDataPoint(const RangeAllType& area) const
  {
    ChromatogramPeak peak_lt {area.getMinRT(), area.getMinIntensity()}, peak_rb {area.getMaxRT(), area.getMaxIntensity()};
    // reference to the current data
    const auto& chrom = getCurrentChrom();
    const Size index = getCurrentIndex();

    // get iterator on first peak with lower position than interval_start
    auto left_it = lower_bound(chrom.begin(), chrom.end(), peak_lt, ChromatogramPeak::PositionLess());

    // get iterator on first peak with higher position than interval_end
    auto right_it = lower_bound(left_it, chrom.end(), peak_rb, ChromatogramPeak::PositionLess());

    if (left_it == right_it) // both are equal => no peak falls into this interval
    {
      return PeakIndex();
    }

    if (left_it == right_it - 1)
    {
      return PeakIndex(index, left_it - chrom.begin());
    }

    auto nearest_it = left_it;
    const auto center_intensity = (peak_lt.getIntensity() + peak_rb.getIntensity()) * 0.5;
    for (auto it = left_it; it != right_it; ++it)
    {
      if (abs(center_intensity - it->getIntensity()) < abs(center_intensity - nearest_it->getIntensity()))
      {
        nearest_it = it;
      }
    }
    return PeakIndex(index, nearest_it - chrom.begin());
  }

  std::unique_ptr<Painter1DBase> LayerData1DChrom::getPainter1D() const
  {
    return make_unique<Painter1DChrom>(this);
  }

  Annotation1DItem* LayerData1DChrom::addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color)
  {
    auto peak = getCurrentChrom()[peak_index.peak];
    auto* item = new Annotation1DPeakItem<decltype(peak)>(peak, text, color);
    item->setSelected(false);
    getCurrentAnnotations().push_front(item);
    return item;
  }

} // namespace OpenMS
