// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/Painter1DBase.h>

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DVerticalLineItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>

// preprocessing and filtering for automated m/z annotations
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>


#include <QPainter>
#include <QPen>

using namespace std;

namespace OpenMS
{
  void Painter1DBase::drawDashedLine(const QPoint& from, const QPoint& to, QPainter* painter, QColor color)
  {
    QPen pen;
    QVector<qreal> dashes;
    dashes << 5 << 5 << 1 << 5;
    pen.setDashPattern(dashes);
    pen.setColor(color);
    painter->save();
    painter->setPen(pen);
    painter->drawLine(from, to);
    painter->restore();
  }

  Painter1DPeak::Painter1DPeak(const LayerDataPeak* parent) : layer_(parent)
  {
  }

  void Painter1DPeak::paint(QPainter* painter, Plot1DCanvas* canvas, int layer_index)
  {
    if (!layer_->visible)
    {
      return;
    }                              
    
    const auto& spectrum = layer_->getCurrentSpectrum();

    // get default icon and peak color
    //QPen icon_pen = QPen(QColor(String(layer_->param.getValue("icon_color").toString()).toQString()), 1);
    QPen pen(QColor(String(layer_->param.getValue("peak_color").toString()).toQString()), 1);
    Qt::PenStyle pen_style = canvas->peak_penstyle_[layer_index];
    pen.setStyle(pen_style);

    painter->setPen(pen);
    // draw dashed elongations for pairs of peaks annotated with a distance
    QColor color = String(canvas->param_.getValue("highlighted_peak_color").toString()).toQString();
    for (auto& it : layer_->getCurrentAnnotations())
    {
      Annotation1DDistanceItem* distance_item = dynamic_cast<Annotation1DDistanceItem*>(it);
      if (distance_item)
      {
        QPoint from, to;
        canvas->dataToWidget(distance_item->getStartPoint().getX(), 0, from, layer_->flipped);

        canvas->dataToWidget(distance_item->getStartPoint().getX(), canvas->visible_area_.maxY(), to, layer_->flipped);
        drawDashedLine(from, to, painter, color);

        canvas->dataToWidget(distance_item->getEndPoint().getX(), 0, from, layer_->flipped);

        canvas->dataToWidget(distance_item->getEndPoint().getX(), canvas->visible_area_.maxY(), to, layer_->flipped);
        drawDashedLine(from, to, painter, color);
      }
    }
   
    auto vbegin = spectrum.MZBegin(canvas->visible_area_.minX());
    auto vend = spectrum.MZEnd(canvas->visible_area_.maxX());
    QPoint begin, end;
    Plot1DCanvas::DrawModes line_mode = canvas->draw_modes_[layer_index];
    switch (line_mode)
    {
      case Plot1DCanvas::DrawModes::DM_PEAKS:
      {
        //---------------------DRAWING PEAKS---------------------

        for (auto it = vbegin; it != vend; ++it)
        {
          if (!layer_->filters.passes(spectrum, it - spectrum.begin())) continue;
          
          // use peak colors stored in the layer, if available
          if (layer_->peak_colors_1d.size() == spectrum.size())
          {
            // find correct peak index
            Size peak_index = std::distance(spectrum.begin(), it);
            pen.setColor(layer_->peak_colors_1d[peak_index]);
            painter->setPen(pen);
          }

          // Warn if non-empty peak color array present but size doesn't match number of peaks
          // This indicates a bug but we gracefully just issue a warning
          if (!layer_->peak_colors_1d.empty() && layer_->peak_colors_1d.size() < spectrum.size())
          {
            OPENMS_LOG_ERROR << "Peak color array size (" << layer_->peak_colors_1d.size() << ") doesn't match number of peaks (" << spectrum.size() << ") in spectrum." << endl;
          }
          canvas->dataToWidget(*it, end, layer_->flipped);
          canvas->dataToWidget(it->getMZ(), 0.0f, begin, layer_->flipped);
          // draw peak
          painter->drawLine(begin, end);
        }
        break;
      }
      case Plot1DCanvas::DrawModes::DM_CONNECTEDLINES:
      {
        //---------------------DRAWING CONNECTED LINES---------------------

        QPainterPath path;

        // connect peaks in visible area; (no clipping needed)
        bool first_point = true;
        for (auto it = vbegin; it != vend; it++)
        {
          canvas->dataToWidget(*it, begin, layer_->flipped);

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

        // clipping on left side
        if (vbegin != spectrum.begin() && vbegin != spectrum.end())
        {
          canvas->dataToWidget(*(vbegin - 1), begin, layer_->flipped);
          canvas->dataToWidget(*(vbegin), end, layer_->flipped);
          painter->drawLine(begin, end);
        }

        // clipping on right side
        if (vend != spectrum.end() && vend != spectrum.begin())
        {
          canvas->dataToWidget(*(vend - 1), begin, layer_->flipped);
          canvas->dataToWidget(*(vend), end, layer_->flipped);
          painter->drawLine(begin, end);
        }

        break;
      }

      default:
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    // annotate interesting m/z's
    if (canvas->draw_interesting_MZs_)
    {
      drawMZAtInterestingPeaks_(*painter, canvas, vbegin, vend);
    }

    // draw all annotation items
    drawAnnotations_(*painter, canvas);

    // draw a legend
    if (canvas->param_.getValue("show_legend").toBool())
    {
      double xpos = canvas->getVisibleArea().maxX() - (canvas->getVisibleArea().maxX() - canvas->getVisibleArea().minX()) * 0.1;
      auto tmp = max_element(spectrum.MZBegin(canvas->visible_area_.minX()), spectrum.MZEnd(xpos), Plot1DCanvas::PeakType::IntensityLess());
      if (tmp != spectrum.end())
      {
        Plot1DCanvas::PointType position(xpos, std::max<double>(tmp->getIntensity() - 100, tmp->getIntensity() * 0.8));
        Annotation1DPeakItem item = Annotation1DPeakItem(position, layer_->getName().toQString(), QColor(String(layer_->param.getValue("peak_color").toString()).toQString()));
        item.draw(canvas, *painter);
      }
    }
  }
  


  void Painter1DPeak::drawAnnotations_(QPainter& painter, Plot1DCanvas* canvas)
  {
    QColor col {QColor(String(layer_->param.getValue("annotation_color").toString()).toQString())};
    // 0: default pen; 1: selected pen
    QPen pen[2] = {col, col.lighter()};

    for (const auto& c : layer_->getCurrentAnnotations())
    {
      painter.setPen(pen[c->isSelected()]);
      c->draw(canvas, painter, layer_->flipped);
    }
  }

  void Painter1DPeak::drawMZAtInterestingPeaks_(QPainter& painter, Plot1DCanvas* canvas, MSSpectrum::ConstIterator vbegin, MSSpectrum::ConstIterator vend)
  {
    if (vbegin == vend)
    {
      return;
    }
    // find interesting peaks

    // copy visible peaks into spec
    MSSpectrum spec;
    for (auto it(vbegin); it != vend; ++it)
    {
      spec.push_back(*it);
    }

    // calculate distance between first and last peak
    --vend;
    double visible_range = vend->getMZ() - vbegin->getMZ();

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

    NLargest nlargest_filter(10); // maximum number of annotated m/z's in visible area
    nlargest_filter.filterPeakSpectrum(spec);
    spec.sortByPosition(); // nlargest changes order

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

      Annotation1DPeakItem item({mz, intensity}, label, Qt::darkGray);
      item.setSelected(false);
      item.draw(canvas, painter, layer_->flipped);
    }
  }

} // namespace OpenMS