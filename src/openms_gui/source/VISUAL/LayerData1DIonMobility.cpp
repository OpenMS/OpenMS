// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerData1DIonMobility.h>

#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/Painter1DBase.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

#include <QMenu>
using namespace std;

namespace OpenMS
{

  std::unique_ptr<LayerStoreData> LayerData1DIonMobility::storeVisibleData(const RangeAllType& /*visible_range*/, const DataFilters& /*layer_filters*/) const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    // does not exist yet...
  }

  std::unique_ptr<LayerStoreData> LayerData1DIonMobility::storeFullData() const
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    // does not exist yet...
  }

  QMenu* LayerData1DIonMobility::getContextMenuAnnotation(Annotation1DItem* /*annot_item*/, bool& /*need_repaint*/)
  {
    auto* context_menu = new QMenu("MobilityPeak", nullptr);
    return context_menu;
  }

  PeakIndex LayerData1DIonMobility::findClosestDataPoint(const RangeAllType& area) const
  {
    MobilityPeak1D peak_lt(area.getMinMobility(), area.getMinIntensity()), peak_rb(area.getMaxMobility(), area.getMaxIntensity());
    // reference to the current data
    const auto& mg = getCurrentMobilogram();
    const Size spectrum_index = getCurrentIndex();

    // get iterator on first peak with lower position than interval_start
    auto left_it = lower_bound(mg.begin(), mg.end(), peak_lt, PeakType::PositionLess());

    // get iterator on first peak with higher position than interval_end
    auto right_it = lower_bound(left_it, mg.end(), peak_rb, PeakType::PositionLess());

    if (left_it == right_it) // both are equal => no peak falls into this interval
    {
      return PeakIndex();
    }

    if (left_it == right_it - 1)
    {
      return PeakIndex(spectrum_index, left_it - mg.begin());
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
    return PeakIndex(spectrum_index, nearest_it - mg.begin());
  }

  std::unique_ptr<Painter1DBase> LayerData1DIonMobility::getPainter1D() const
  {
    return make_unique<Painter1DIonMobility>(this);
  }

  Annotation1DItem* LayerData1DIonMobility::addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color)
  {
    auto peak = getCurrentMobilogram()[peak_index.peak];
    auto* item = new Annotation1DPeakItem<decltype(peak)>(peak, text, color);
    item->setSelected(false);
    getCurrentAnnotations().push_front(item);
    return item;
  }

  


}// namespace OpenMS
