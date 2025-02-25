// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerData1DChrom.h>
#include <OpenMS/VISUAL/LayerData1DIonMobility.h>
#include <OpenMS/VISUAL/LayerData1DPeak.h>
#include <OpenMS/VISUAL/Painter1DBase.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>

// preprocessing and filtering for automated m/z annotations
#include <OpenMS/PROCESSING/DEISOTOPING/Deisotoper.h>
#include <OpenMS/PROCESSING/FILTERING/ThresholdMower.h>
#include <OpenMS/PROCESSING/FILTERING/NLargest.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>


#include <QPainter>
#include <QPen>

using namespace std;

namespace OpenMS
{
  void Painter1DBase::drawAnnotations_(const LayerData1DBase* layer, QPainter& painter, Plot1DCanvas* canvas) const
  {
    const QColor col {QColor(String(layer->param.getValue("annotation_color").toString()).toQString())};
    // 0: default pen; 1: selected pen
    const QPen pen[2] = {col, col.lighter()};

    for (const auto& c : layer->getCurrentAnnotations())
    {
      painter.setPen(pen[c->isSelected()]);
      c->draw(canvas, painter, layer->flipped);
    }
  }

  // ###############################
  // ###### 1D Peak
  // ###############################

  Painter1DPeak::Painter1DPeak(const LayerData1DPeak* parent) : layer_(parent)
  {
  }

  void Painter1DPeak::paint(QPainter* painter, Plot1DCanvas* canvas, int layer_index)
  {
    if (!layer_->visible)
    {
      return;
    }                              
    
    const auto& spectrum = layer_->getCurrentSpectrum();

    // get default peak color
    QPen pen(QColor(String(layer_->param.getValue("peak_color").toString()).toQString()), 1);
    pen.setStyle(canvas->peak_penstyle_[layer_index]);
    painter->setPen(pen);

    // draw dashed elongations for pairs of peaks annotated with a distance
    const QColor color = String(canvas->param_.getValue("highlighted_peak_color").toString()).toQString();
    for (auto& it : layer_->getCurrentAnnotations())
    {
      const auto distance_item = dynamic_cast<Annotation1DDistanceItem*>(it);
      if (!distance_item) continue;

      auto draw_line_ = [&](const PointXYType& p) {
        QPoint from;
        canvas->dataToWidget(p, from, layer_->flipped);
        from = canvas->getGravitator().gravitateZero(from);
        QPoint to = canvas->getGravitator().gravitateMax(from, canvas->canvasPixelArea());
        drawDashedLine(from, to, painter, color);
      };
      draw_line_(distance_item->getStartPoint());
      draw_line_(distance_item->getEndPoint());
    }
   
    const auto v_begin = spectrum.MZBegin(canvas->visible_area_.getAreaUnit().getMinMZ());
    const auto v_end = spectrum.MZEnd(canvas->visible_area_.getAreaUnit().getMaxMZ());
    QPoint begin, end;
    
    switch (canvas->draw_modes_[layer_index])
    {
      case Plot1DCanvas::DrawModes::DM_PEAKS:
      {
        //---------------------DRAWING PEAKS---------------------

        for (auto it = v_begin; it != v_end; ++it)
        {
          if (!layer_->filters.passes(spectrum, it - spectrum.begin())) continue;
          
          // use peak colors stored in the layer, if available
          if (layer_->peak_colors_1d.size() == spectrum.size())
          {
            // find correct peak index
            const Size peak_index = std::distance(spectrum.cbegin(), it);
            pen.setColor(layer_->peak_colors_1d[peak_index]);
            painter->setPen(pen);
          }
          else if (!layer_->peak_colors_1d.empty())
          { // Warn if non-empty peak color array present but size doesn't match number of peaks
            // This indicates a bug but we gracefully just issue a warning
            OPENMS_LOG_ERROR << "Peak color array size (" << layer_->peak_colors_1d.size() << ") doesn't match number of peaks (" << spectrum.size() << ") in spectrum." << endl;
          }
          // draw stick
          auto p_xy = canvas->getMapper().map(*it);
          canvas->dataToWidget(p_xy, end, layer_->flipped);
          canvas->dataToWidget(canvas->getGravitator().gravitateZero(p_xy), begin, layer_->flipped);
          painter->drawLine(begin, end);
        }
        break;
      }
      case Plot1DCanvas::DrawModes::DM_CONNECTEDLINES:
      {
        //---------------------DRAWING CONNECTED LINES---------------------

        QPainterPath path;

        // connect peaks in visible area; 
        // clipping on left and right side
        auto v_begin_cl = v_begin;
        if (v_begin_cl != spectrum.cbegin() && v_begin_cl != spectrum.cend())
          --v_begin_cl;
        auto v_end_cl = v_end;
        if (v_end_cl != spectrum.cbegin() && v_end_cl != spectrum.cend())
          ++v_end_cl;

        bool first_point = true;
        for (auto it = v_begin_cl; it != v_end_cl; ++it)
        {
          if (!layer_->filters.passes(spectrum, it - spectrum.begin())) continue;

          canvas->dataToWidget(canvas->getMapper().map(*it), begin, layer_->flipped);

          // connect lines
          if (first_point)
          {
            path.moveTo(begin);
            first_point = false;
          }
          else
          {
            path.lineTo(begin);
          }
        }
        painter->drawPath(path);
        break;
      }

      default:
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    // annotate interesting m/z's
    if (canvas->draw_interesting_MZs_)
    {
      drawMZAtInterestingPeaks_(*painter, canvas, v_begin, v_end);
    }

    // draw all annotation items
    drawAnnotations_(layer_, *painter, canvas);
  }

  void Painter1DPeak::drawMZAtInterestingPeaks_(QPainter& painter, Plot1DCanvas* canvas, MSSpectrum::ConstIterator v_begin, MSSpectrum::ConstIterator v_end) const
  {
    if (v_begin == v_end)
    {
      return;
    }
    // find interesting peaks

    // copy visible peaks into spec
    MSSpectrum spec;
    for (auto it(v_begin); it != v_end; ++it)
    {
      spec.push_back(*it);
    }

    // calculate distance between first and last peak
    --v_end;
    double visible_range = v_end->getMZ() - v_begin->getMZ();

    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakSpectrum(spec);

    // deisotope so we don't consider higher isotopic peaks
    Deisotoper::deisotopeAndSingleCharge(spec,
                                         100,   // tolerance
                                         true,  // ppm
                                         1, 6,  // min / max charge
                                         false, // keep only deisotoped
                                         3, 10, // min / max isopeaks
                                         false, // don't convert fragment m/z to mono-charge
                                         true); // annotate charge in integer data array

    // filter for local high-intensity peaks
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    double window_size = visible_range / 10.0;
    filter_param.setValue("windowsize", window_size, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 2, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "slide", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);
    window_mower_filter.filterPeakSpectrum(spec);

    // maximum number of annotated m/z's in visible area
    NLargest(10).filterPeakSpectrum(spec);
    spec.sortByPosition(); // n-largest changes order

    for (size_t i = 0; i < spec.size(); ++i)
    {
      auto mz(spec[i].getMZ());
      auto intensity(spec[i].getIntensity());

      QString label = String::number(mz, 4).toQString();

      if (!spec.getIntegerDataArrays().empty() && spec.getIntegerDataArrays()[0].size() == spec.size())
      {
        int charge = spec.getIntegerDataArrays()[0][i];
        // TODO: handle negative mode

        // here we explicitly also annotate singly charged ions to distinguish them from unknown charge (0)
        if (charge != 0)
        {
          label += charge == 1 ? "<sup>+</sup>" : "<sup>" + QString::number(charge) + "+</sup>";
        }
      }

      Annotation1DPeakItem item(Peak1D{mz, intensity}, label, Qt::darkGray);
      item.setSelected(false);
      item.draw(canvas, painter, layer_->flipped);
    }
  }


  // ###############################
  // ###### 1D Chrom
  // ###############################

  Painter1DChrom::Painter1DChrom(const LayerData1DChrom* parent) : layer_(parent)
  {
  }

  void Painter1DChrom::paint(QPainter* painter, Plot1DCanvas* canvas, int layer_index)
  {
    if (!layer_->visible)
    {
      return;
    }

    const auto& data = layer_->getCurrentChrom();

    // get default peak color
    QPen pen(QColor(String(layer_->param.getValue("peak_color").toString()).toQString()), 1);
    pen.setStyle(canvas->peak_penstyle_[layer_index]);
    painter->setPen(pen);

    const auto v_begin = data.RTBegin(canvas->visible_area_.getAreaUnit().getMinRT());
    const auto v_end = data.RTEnd(canvas->visible_area_.getAreaUnit().getMaxRT());
    QPoint begin, end;
    switch (canvas->draw_modes_[layer_index])
    {
      case Plot1DCanvas::DrawModes::DM_PEAKS: {
        //---------------------DRAWING PEAKS---------------------

        for (auto it = v_begin; it != v_end; ++it)
        {
          if (!layer_->filters.passes(data, it - data.begin()))
            continue;

          // use peak colors stored in the layer, if available
          if (layer_->peak_colors_1d.size() == data.size())
          {
            // find correct peak index
            const Size peak_index = std::distance(data.begin(), it);
            pen.setColor(layer_->peak_colors_1d[peak_index]);
            painter->setPen(pen);
          }
          else if (!layer_->peak_colors_1d.empty())
          { // Warn if non-empty peak color array present but size doesn't match number of peaks
            // This indicates a bug but we gracefully just issue a warning
            OPENMS_LOG_ERROR << "Peak color array size (" << layer_->peak_colors_1d.size() << ") doesn't match number of peaks (" << data.size() << ") in chromatogram." << endl;
          }
          // draw stick
          auto p_xy = canvas->getMapper().map(*it);
          canvas->dataToWidget(p_xy, end, layer_->flipped);
          canvas->dataToWidget(canvas->getGravitator().gravitateZero(p_xy), begin, layer_->flipped);
          painter->drawLine(begin, end);
        }
        break;
      }
      case Plot1DCanvas::DrawModes::DM_CONNECTEDLINES: {
        //---------------------DRAWING CONNECTED LINES---------------------

        QPainterPath path;

        // connect peaks in visible area;
        // clipping on left and right side
        auto v_begin_cl = v_begin;
        if (v_begin_cl != data.cbegin() && v_begin_cl != data.cend())
          --v_begin_cl;
        auto v_end_cl = v_end;
        if (v_end_cl != data.cbegin() && v_end_cl != data.cend())
          ++v_end_cl;

        bool first_point = true;
        for (auto it = v_begin_cl; it != v_end_cl; ++it)
        {
          if (!layer_->filters.passes(data, it - data.begin())) continue;

          canvas->dataToWidget(canvas->getMapper().map(*it), begin, layer_->flipped);

          // connect lines
          if (first_point)
          {
            path.moveTo(begin);
            first_point = false;
          }
          else
          {
            path.lineTo(begin);
          }
        }
        painter->drawPath(path);
        break;
      }

      default:
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    // draw all annotation items
    drawAnnotations_(layer_, *painter, canvas);
  }

  // ###############################
  // ###### 1D Mobilogram
  // ###############################

  Painter1DIonMobility::Painter1DIonMobility(const LayerData1DIonMobility* parent) : layer_(parent)
  {
  }

  void Painter1DIonMobility::paint(QPainter* painter, Plot1DCanvas* canvas, int layer_index)
  {
    if (!layer_->visible)
    {
      return;
    }

    const auto& data = layer_->getCurrentMobilogram();

    // get default peak color
    QPen pen(QColor(String(layer_->param.getValue("peak_color").toString()).toQString()), 1);
    pen.setStyle(canvas->peak_penstyle_[layer_index]);
    painter->setPen(pen);

    const auto v_begin = data.MBBegin(canvas->visible_area_.getAreaUnit().getMinMobility());
    const auto v_end = data.MBEnd(canvas->visible_area_.getAreaUnit().getMaxMobility());
    QPoint begin, end;
    switch (canvas->draw_modes_[layer_index])
    {
      case Plot1DCanvas::DrawModes::DM_PEAKS: {
        //---------------------DRAWING PEAKS---------------------

        for (auto it = v_begin; it != v_end; ++it)
        {
          //if (!layer_->filters.passes(data, it - data.begin()))
          //  continue;

          // use peak colors stored in the layer, if available
          if (layer_->peak_colors_1d.size() == data.size())
          {
            // find correct peak index
            const Size peak_index = std::distance(data.begin(), it);
            pen.setColor(layer_->peak_colors_1d[peak_index]);
            painter->setPen(pen);
          }
          else if (!layer_->peak_colors_1d.empty())
          { // Warn if non-empty peak color array present but size doesn't match number of peaks
            // This indicates a bug but we gracefully just issue a warning
            OPENMS_LOG_ERROR << "Peak color array size (" << layer_->peak_colors_1d.size() << ") doesn't match number of peaks (" << data.size() << ") in chromatogram." << endl;
          }
          // draw stick
          auto p_xy = canvas->getMapper().map(*it);
          canvas->dataToWidget(p_xy, end, layer_->flipped);
          canvas->dataToWidget(canvas->getGravitator().gravitateZero(p_xy), begin, layer_->flipped);
          painter->drawLine(begin, end);
        }
        break;
      }
      case Plot1DCanvas::DrawModes::DM_CONNECTEDLINES: {
        //---------------------DRAWING CONNECTED LINES---------------------

        QPainterPath path;

        // connect peaks in visible area;
        // clipping on left and right side
        auto v_begin_cl = v_begin;
        if (v_begin_cl != data.cbegin() && v_begin_cl != data.cend())
          --v_begin_cl;
        auto v_end_cl = v_end;
        if (v_end_cl != data.cbegin() && v_end_cl != data.cend())
          ++v_end_cl;

        bool first_point = true;
        for (auto it = v_begin_cl; it != v_end_cl; ++it)
        {
          if (!layer_->filters.passes(data, it - data.begin())) continue;
          
          canvas->dataToWidget(canvas->getMapper().map(*it), begin, layer_->flipped);

          // connect lines
          if (first_point)
          {
            path.moveTo(begin);
            first_point = false;
          }
          else
          {
            path.lineTo(begin);
          }
        }
        painter->drawPath(path);
        break;
      }

      default:
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    // draw all annotation items
    drawAnnotations_(layer_, *painter, canvas);
  }
} // namespace OpenMS