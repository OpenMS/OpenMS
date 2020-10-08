// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

// Qt
#include <QMouseEvent>
#include <QtWidgets/QMessageBox>
#include <QPainterPath>
#include <QPainter>
#include <QtCore/QTime>
#include <QtWidgets/QMenu>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QInputDialog>
#include <QtSvg/QSvgGenerator>

// OpenMS
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DVerticalLineItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DPrefDialog.h>

// preprocessing and filtering for automated m/z annotations
#include <OpenMS/FILTERING/DATAREDUCTION/Deisotoper.h>
#include <OpenMS/FILTERING/TRANSFORMERS/ThresholdMower.h>
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>
#include <OpenMS/FILTERING/TRANSFORMERS/WindowMower.h>

#include <iostream>
#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

namespace OpenMS
{
  using namespace Math;
  using namespace Internal;

  Spectrum1DCanvas::Spectrum1DCanvas(const Param& preferences, QWidget* parent) :
    SpectrumCanvas(preferences, parent),
    mirror_mode_(false),
    moving_annotations_(false),
    show_alignment_(false),
    aligned_peaks_mz_delta_(),
    alignment_score_(0),
    is_swapped_(true),
    ion_ladder_visible_(true),
    draw_interesting_MZs_(false)
  {
    //Parameter handling
    defaults_.setValue("highlighted_peak_color", "#ff0000", "Highlighted peak color.");
    defaults_.setValue("icon_color", "#000000", "Peak icon color.");
    defaults_.setValue("peak_color", "#0000ff", "Peak color.");
    defaults_.setValue("annotation_color", "#000055", "Annotation color.");
    defaults_.setValue("background_color", "#ffffff", "Background color.");
    defaults_.setValue("show_legend", "false", "Annotate each layer with its name on the canvas.");
    defaultsToParam_();
    setName("Spectrum1DCanvas");
    setParameters(preferences);

    //connect preferences change to the right slot
    connect(this, SIGNAL(preferencesChange()), this, SLOT(currentLayerParamtersChanged_()));
  }

  Spectrum1DCanvas::~Spectrum1DCanvas()
  {
  }

  bool Spectrum1DCanvas::addChromLayer(ExperimentSharedPtrType chrom_exp_sptr, ODExperimentSharedPtrType ondisc_sptr, const String& filename,
                                       const String& caption,
                                       ExperimentSharedPtrType exp_sptr,
                                       const int index,
                                       const bool multiple_select)
  {
    // we do not want addLayer to trigger repaint, since we have not set the chromatogram data!
    this->blockSignals(true);
    RAIICleanup clean([&]()
    {
      this->blockSignals(false);
    });

    // add chromatogram data as peak spectrum
    if (!addLayer(chrom_exp_sptr, ondisc_sptr, filename))
    {
      return false;
    }

    setDrawMode(Spectrum1DCanvas::DM_CONNECTEDLINES);
    setIntensityMode(Spectrum1DCanvas::IM_NONE);

    getCurrentLayer().setName(caption);
    getCurrentLayer().getChromatogramData() = exp_sptr; // save the original chromatogram data so that we can access it later
    //this is a hack to store that we have chromatogram data, that we selected multiple ones and which one we selected
    getCurrentLayer().getPeakDataMuteable()->setMetaValue("is_chromatogram", "true");
    getCurrentLayer().getPeakDataMuteable()->setMetaValue("multiple_select", multiple_select ? "true" : "false");
    getCurrentLayer().getPeakDataMuteable()->setMetaValue("selected_chromatogram", index);

    return true;
  }

  void Spectrum1DCanvas::activateLayer(Size layer_index)
  {
    layers_.setCurrentLayer(layer_index);

    // no peak is selected
    selected_peak_.clear();

    emit layerActivated(this);
  }

  void Spectrum1DCanvas::setVisibleArea(DRange<2> range)
  {
    changeVisibleArea_(AreaType(range.minX(), visible_area_.minY(), range.maxX(), visible_area_.maxY()));
  }

  void Spectrum1DCanvas::changeVisibleArea_(double lo, double hi, bool repaint, bool add_to_stack)
  {
    changeVisibleArea_(AreaType(lo, visible_area_.minY(), hi, visible_area_.maxY()), repaint, add_to_stack);
    emit layerZoomChanged(this);
  }

  void Spectrum1DCanvas::dataToWidget(const PeakType& peak, QPoint& point, bool flipped, bool percentage)
  {
    dataToWidget(peak.getMZ(), peak.getIntensity(), point, flipped, percentage);
  }

  void Spectrum1DCanvas::dataToWidget(double x, double y, QPoint& point, bool flipped, bool percentage)
  {
    QPoint tmp;
    if (percentage)
    {
      y *= getSnapFactor() * percentage_factor_;
    }
    SpectrumCanvas::dataToWidget_(x, y, tmp);
    point.setX(tmp.x());
    double alignment_shrink_factor = 1.0;
    if (height() > 10)
    {
      alignment_shrink_factor = (double)(height() - 10) / (double)height();
    }
    if (mirror_mode_)
    {
      if (flipped)
      {
        if (!show_alignment_)
        {
          point.setY(height() - (int)(tmp.y() / 2.0));
        }
        else         // show_alignment_
        {
          point.setY(height() - (int)((tmp.y() * alignment_shrink_factor) / 2.0));
        }
      }
      else       // !flipped
      {
        if (!show_alignment_)
        {
          point.setY((int)(tmp.y() / 2.0));
        }
        else         // show_alignment_
        {
          point.setY((int)((tmp.y() * alignment_shrink_factor) / 2.0));
        }
      }
    }
    else     // !mirror_mode_
    {
      point.setY((int)(tmp.y()));
    }
  }

  SpectrumCanvas::PointType Spectrum1DCanvas::widgetToData(const QPoint& pos, bool percentage)
  {
    return widgetToData(pos.x(), pos.y(), percentage);
  }

  SpectrumCanvas::PointType Spectrum1DCanvas::widgetToData(double x, double y, bool percentage)
  {
    double actual_y;
    double alignment_shrink_factor = 1.0;
    if (height() > 10)
    {
      alignment_shrink_factor = (double)(height() - 10) / (double)height();
    }

    if (mirror_mode_)
    {
      if (y > height() / 2.0)
      {
        if (!show_alignment_)
        {
          actual_y = (height() - y) * 2;
        }
        else
        {
          actual_y = (height() - y) * 2 / alignment_shrink_factor;
        }
      }
      else       // y <= height()/2
      {
        if (!show_alignment_)
        {
          actual_y = y * 2;
        }
        else
        {
          actual_y = y * 2 / alignment_shrink_factor;
        }
      }
    }
    else
    {
      actual_y = y;
    }
    PointType p = SpectrumCanvas::widgetToData_(x, actual_y);
    if (percentage)
    {
      p.setY(p.getY() / (getSnapFactor() * percentage_factor_));
    }
    return p;
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Qt events

  void Spectrum1DCanvas::mousePressEvent(QMouseEvent* e)
  {
    // get mouse position in widget coordinates
    last_mouse_pos_ = e->pos();

    if (e->button() == Qt::LeftButton)
    {
      // selection/deselection of annotation items
      Annotation1DItem* item = getCurrentLayer().getCurrentAnnotations().getItemAt(last_mouse_pos_);
      if (item)
      {
        if (!(e->modifiers() & Qt::ControlModifier))
        {
          // edit via double-click
          if (e->type() == QEvent::MouseButtonDblClick)
          {
            item->editText();
          }
          else if (!item->isSelected())
          {
            // the item becomes the only selected item
            getCurrentLayer().getCurrentAnnotations().deselectAll();
            item->setSelected(true);
          }
          // an item was clicked -> can be moved on the canvas
          moving_annotations_ = true;
        }
        else
        {
          // ctrl pressed -> allow selection/deselection of multiple items, do not deselect others
          item->setSelected(!item->isSelected());
        }

        // if item is a distance item: show distance of selected item in status bar
        Annotation1DDistanceItem * distance_item = dynamic_cast<Annotation1DDistanceItem *>(item);
        if (distance_item)
        {
          const double start_p = distance_item->getStartPoint().getX();
          const double end_p = distance_item->getEndPoint().getX();
          emit sendStatusMessage(QString("Measured: dMZ = %1").arg(end_p - start_p).toStdString(), 0);
        }
      }
      else
      {
        // no item was under the cursor
        getCurrentLayer().getCurrentAnnotations().deselectAll();
      }

      if (action_mode_ == AM_ZOOM)
      {
        rubber_band_.setGeometry(QRect(e->pos(), QSize()));
        rubber_band_.show();
      }
      else if (action_mode_ == AM_MEASURE)
      {
        if (isMzToXAxis())
        {
          if (selected_peak_.isValid())
          {
            measurement_start_ = selected_peak_;
            const ExperimentType::PeakType & peak = getCurrentLayer().getCurrentSpectrum()[measurement_start_.peak];
            if (intensity_mode_ == IM_PERCENTAGE)
            {
              updatePercentageFactor_(getCurrentLayerIndex());
            }
            else
            {
              percentage_factor_ = 1.0;
            }
            dataToWidget(peak, measurement_start_point_, getCurrentLayer().flipped);
            measurement_start_point_.setY(last_mouse_pos_.y());
          }
          else
          {
            measurement_start_.clear();
          }
        }
        else         // !isMzToXAxis()
        {
          if (selected_peak_.isValid())
          {
            measurement_start_ = selected_peak_;
            const ExperimentType::PeakType & peak = getCurrentLayer().getCurrentSpectrum()[measurement_start_.peak];
            updatePercentageFactor_(getCurrentLayerIndex());
            dataToWidget(peak, measurement_start_point_, getCurrentLayer().flipped);
            measurement_start_point_.setX(last_mouse_pos_.x());
          }
          else
          {
            measurement_start_.clear();
          }
        }
      }
    }
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum1DCanvas::mouseMoveEvent(QMouseEvent* e)
  {
    // mouse position relative to the diagram widget
    QPoint p = e->pos();
    PointType data_pos = widgetToData(p);
    emit sendCursorStatus(data_pos.getX(), getCurrentLayer().getCurrentSpectrum().getRT());

    PeakIndex near_peak = findPeakAtPosition_(p);

    if (e->buttons() & Qt::LeftButton)
    {
      bool move = moving_annotations_;
      if (mirror_mode_ && (getCurrentLayer().flipped ^ (p.y() > height() / 2)))
      {
        move = false;
      }
      if (move)
      {
        updatePercentageFactor_(getCurrentLayerIndex());
        PointType delta = widgetToData(p, true) - widgetToData(last_mouse_pos_, true);

        Annotations1DContainer& ann_1d = getCurrentLayer().getCurrentAnnotations();
        for (Annotations1DContainer::Iterator it = ann_1d.begin(); it != ann_1d.end(); ++it)
        {
          if ((*it)->isSelected())
          {
            (*it)->move(delta);
          }
        }
        update_(OPENMS_PRETTY_FUNCTION);
        last_mouse_pos_ = p;
      }
      else if (action_mode_ == AM_TRANSLATE)
      {
        // translation in data metric
        double shift = widgetToData(last_mouse_pos_).getX() - widgetToData(p).getX();
        double newLo = visible_area_.minX() + shift;
        double newHi = visible_area_.maxX() + shift;
        // check if we are falling out of bounds
        if (newLo < overall_data_range_.minX())
        {
          newLo = overall_data_range_.minX();
          newHi = newLo + visible_area_.width();
        }
        if (newHi > overall_data_range_.maxX())
        {
          newHi = overall_data_range_.maxX();
          newLo = newHi - visible_area_.width();
        }
        //chage data area
        changeVisibleArea_(newLo, newHi);
        last_mouse_pos_ = p;
      }
      else if (action_mode_ == AM_MEASURE)
      {
        if (near_peak.peak != measurement_start_.peak)
        {
          selected_peak_ = near_peak;
          last_mouse_pos_ = p;
          update_(OPENMS_PRETTY_FUNCTION);
        }
      }
      else if (action_mode_ == AM_ZOOM)
      {
        PointType pos = widgetToData(p);

        if (isMzToXAxis())
        {
          rubber_band_.setGeometry(QRect(last_mouse_pos_.x(), 0, p.x() - last_mouse_pos_.x(), height()).normalized());
        }
        else
        {
          rubber_band_.setGeometry(QRect(0, last_mouse_pos_.y(), width(), p.y() - last_mouse_pos_.y()).normalized());
        }
        rubber_band_.show();         //if the mouse button is pressed before the zoom key is pressed

        update_(OPENMS_PRETTY_FUNCTION);
      }
    }
    else if (!e->buttons())     //no buttons pressed
    {
      selected_peak_ = findPeakAtPosition_(p);
      update_(OPENMS_PRETTY_FUNCTION);
    }

    //show coordinates
    if (selected_peak_.isValid())
    {
      String status;
      const ExperimentType::SpectrumType& s = getCurrentLayer().getCurrentSpectrum();
      for (Size m = 0; m < s.getFloatDataArrays().size(); ++m)
      {
        if (selected_peak_.peak < s.getFloatDataArrays()[m].size())
        {
          status += s.getFloatDataArrays()[m].getName() + ": " + s.getFloatDataArrays()[m][selected_peak_.peak] + " ";
        }
      }
      for (Size m = 0; m < s.getIntegerDataArrays().size(); ++m)
      {
        if (selected_peak_.peak < s.getIntegerDataArrays()[m].size())
        {
          status += s.getIntegerDataArrays()[m].getName() + ": " + s.getIntegerDataArrays()[m][selected_peak_.peak] + " ";
        }
      }
      for (Size m = 0; m < s.getStringDataArrays().size(); ++m)
      {
        if (selected_peak_.peak < s.getStringDataArrays()[m].size())
        {
          status += s.getStringDataArrays()[m].getName() + ": " + s.getStringDataArrays()[m][selected_peak_.peak] + " ";
        }
      }
      emit sendStatusMessage(status, 0);
    }
  }

  void Spectrum1DCanvas::mouseReleaseEvent(QMouseEvent* e)
  {
    if (e->button() == Qt::LeftButton)
    {
      if (action_mode_ == AM_ZOOM)
      {
        rubber_band_.hide();
        QRect rect = rubber_band_.geometry();
        if (rect.width() != 0)
        {
          AreaType area(widgetToData(rect.topLeft()), widgetToData(rect.bottomRight()));
          changeVisibleArea_(area.minX(), area.maxX(), true, true);
        }
      }
      else if (action_mode_ == AM_MEASURE)
      {
        if (!selected_peak_.isValid())
        {
          measurement_start_.clear();
        }
        if (measurement_start_.isValid() && selected_peak_.peak != measurement_start_.peak)
        {
          const ExperimentType::PeakType& peak_1 = getCurrentLayer().getCurrentSpectrum()[measurement_start_.peak];
          const ExperimentType::PeakType& peak_2 = getCurrentLayer().getCurrentSpectrum()[selected_peak_.peak];
          updatePercentageFactor_(getCurrentLayerIndex());
          PointType p = widgetToData(measurement_start_point_, true);
          bool peak_1_less = peak_1.getMZ() < peak_2.getMZ();
          double start_mz = peak_1_less ? peak_1.getMZ() : peak_2.getMZ();
          double end_mz = peak_1_less ? peak_2.getMZ() : peak_1.getMZ();
          double distance = end_mz - start_mz;
          PointType start_p(start_mz, p.getY());
          PointType end_p(end_mz, p.getY());
          // draw line for measured distance between two peaks and annotate with distance in m/z -- use 4 digits to resolve 13C distances between isotopes
          Annotation1DItem * item = new Annotation1DDistanceItem(QString::number(distance, 'f', 4), start_p, end_p);
          getCurrentLayer().getCurrentAnnotations().push_front(item);
        }
      }

      ensureAnnotationsWithinDataRange_();
      moving_annotations_ = false;

      measurement_start_.clear();
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  void Spectrum1DCanvas::keyPressEvent(QKeyEvent* e)
  {
    // Delete pressed => delete selected annotations from the current layer
    if (e->key() == Qt::Key_Delete)
    {
      e->accept();
      getCurrentLayer().removePeakAnnotationsFromPeptideHit(getCurrentLayer().getCurrentAnnotations().getSelectedItems());
      getCurrentLayer().getCurrentAnnotations().removeSelectedItems();
      update_(OPENMS_PRETTY_FUNCTION);
    }
    // 'a' pressed && in zoom mode (ctrl pressed) => select all annotation items
    else if ((e->modifiers() & Qt::ControlModifier) && (e->key() == Qt::Key_A))
    {
      e->accept();
      getCurrentLayer().getCurrentAnnotations().selectAll();
      update_(OPENMS_PRETTY_FUNCTION);
    }
    else
    {
      SpectrumCanvas::keyPressEvent(e);
    }
  }

  PeakIndex Spectrum1DCanvas::findPeakAtPosition_(QPoint p)
  {
    //no layers => return invalid peak index
    if (layers_.empty())
      return PeakIndex();

    // mirror mode and p not on same half as active layer => return invalid peak index
    if (mirror_mode_ && (getCurrentLayer().flipped ^ (p.y() > height() / 2)))
      return PeakIndex();

    //reference to the current data
    const SpectrumType& spectrum = getCurrentLayer().getCurrentSpectrum();
    Size spectrum_index = getCurrentLayer().getCurrentSpectrumIndex();

    // get the interval (in diagramm metric) that will be projected on screen coordinate p.x() or p.y() (depending on orientation)
    PointType lt = widgetToData(p - QPoint(2, 2), true);
    PointType rb = widgetToData(p + QPoint(2, 2), true);

    // get iterator on first peak with higher position than interval_start
    PeakType temp;
    temp.setMZ(min(lt.getX(), rb.getX()));
    SpectrumConstIteratorType left_it = lower_bound(spectrum.begin(), spectrum.end(), temp, PeakType::PositionLess());

    // get iterator on first peak with higher position than interval_end
    temp.setMZ(max(lt.getX(), rb.getX()));
    SpectrumConstIteratorType right_it = lower_bound(left_it, spectrum.end(), temp, PeakType::PositionLess());

    if (left_it == right_it)     // both are equal => no peak falls into this interval
    {
      return PeakIndex();
    }

    if (left_it == right_it - 1)
    {
      return PeakIndex(spectrum_index, left_it - spectrum.begin());
    }

    SpectrumConstIteratorType nearest_it = left_it;

    // select source interval start and end depending on diagram orientation
    updatePercentageFactor_(getCurrentLayerIndex());
    QPoint tmp;
    dataToWidget(0, overall_data_range_.minY(), tmp, getCurrentLayer().flipped, true);
    double dest_interval_start = tmp.y();
    dataToWidget(0, overall_data_range_.maxY(), tmp, getCurrentLayer().flipped, true);
    double dest_interval_end = tmp.y();

    int nearest_intensity = static_cast<int>(intervalTransformation(nearest_it->getIntensity(), visible_area_.minY(),
                                                                    visible_area_.maxY(), dest_interval_start, dest_interval_end));
    for (SpectrumConstIteratorType it = left_it; it != right_it; it++)
    {
      int current_intensity = static_cast<int>(intervalTransformation(it->getIntensity(), visible_area_.minY(), visible_area_.maxY(),
                                                                  dest_interval_start, dest_interval_end));
      if (abs(current_intensity - p.y()) < abs(nearest_intensity - p.y()))
      {
        nearest_intensity = current_intensity;
        nearest_it = it;
      }
    }
    return PeakIndex(spectrum_index, nearest_it - spectrum.begin());
  }

  //////////////////////////////////////////////////////////////////////////////////
  // SLOTS

  void Spectrum1DCanvas::removeLayer(Size layer_index)
  {
    //remove settings
    layers_.removeLayer(layer_index);
    draw_modes_.erase(draw_modes_.begin() + layer_index);
    peak_penstyle_.erase(peak_penstyle_.begin() + layer_index);

    //update nearest peak
    selected_peak_.clear();

    //abort if there are no layers anymore
    if (layers_.empty())
    {
      overall_data_range_ = DRange<3>::empty;
      update_(OPENMS_PRETTY_FUNCTION);
      return;
    }

    if (!flippedLayersExist())
    {
      setMirrorModeActive(false);
    }

    //update range area
    recalculateRanges_(0, 2, 1);
    double width = overall_data_range_.width();
    overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
    overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
    overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());

    zoomClear_();

    if (overall_data_range_.maxX() - overall_data_range_.minX() < 1.0)
    {
      AreaType new_area(overall_data_range_.minX() - 1.0, overall_data_range_.minY(),
                        overall_data_range_.maxX() + 1.0, overall_data_range_.maxY());
      changeVisibleArea_(new_area, true, true);
    }
    else
    {
      AreaType new_area(overall_data_range_.minX(), overall_data_range_.minY(),
                        overall_data_range_.maxX(), overall_data_range_.maxY());
      changeVisibleArea_(new_area, true, true);
    }
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum1DCanvas::setDrawMode(DrawModes mode)
  {
    //no layers
    if (layers_.empty())
      return;

    if (draw_modes_[getCurrentLayerIndex()] != mode)
    {
      draw_modes_[getCurrentLayerIndex()] = mode;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  Spectrum1DCanvas::DrawModes Spectrum1DCanvas::getDrawMode() const
  {
    //no layers
    if (layers_.empty())
      return DM_PEAKS;

    return draw_modes_[getCurrentLayerIndex()];
  }

  void Spectrum1DCanvas::paintEvent(QPaintEvent* e)
  {
    QPainter painter(this);
    paint(&painter, e);
    painter.end();
  }

  void Spectrum1DCanvas::paint(QPainter* painter, QPaintEvent* e)
  {
    QTime timer;
    timer.start();

    // clear
    painter->fillRect(0, 0, this->width(), this->height(),
                      QColor(param_.getValue("background_color").toQString()));

    // only fill background if no layer is present
    if (getLayerCount() == 0)
    {
      e->accept();
      return;
    }

    // gridlines
    emit recalculateAxes();
    paintGridLines_(*painter);

    for (Size i = 0; i < getLayerCount(); ++i)
    {
      const LayerData& layer = getLayer(i);

      // skip non peak data layer or invisible
      if (layer.type != LayerData::DT_PEAK || !layer.visible) { continue; }

      const ExperimentType::SpectrumType& spectrum = layer.getCurrentSpectrum();

      // get default icon and peak color
      QPen icon_pen = QPen(QColor(layer.param.getValue("icon_color").toQString()), 1);
      QPen pen(QColor(layer.param.getValue("peak_color").toQString()), 1);
      pen.setStyle(peak_penstyle_[i]);

      // TODO option for variable pen width
      // pen.setWidthF(1.5);
      painter->setPen(pen);
      updatePercentageFactor_(i);
      SpectrumConstIteratorType vbegin = getLayer(i).getCurrentSpectrum().MZBegin(visible_area_.minX());
      SpectrumConstIteratorType vend = getLayer(i).getCurrentSpectrum().MZEnd(visible_area_.maxX());

      // draw dashed elongations for pairs of peaks annotated with a distance
      for (auto it = layer.getCurrentAnnotations().begin();
            it != layer.getCurrentAnnotations().end(); ++it)
      {
        Annotation1DDistanceItem* distance_item = dynamic_cast<Annotation1DDistanceItem*>(*it);
        if (distance_item)
        {
          QPoint from, to;
          dataToWidget(distance_item->getStartPoint().getX(), 0, from, layer.flipped);

          dataToWidget(distance_item->getStartPoint().getX(), getVisibleArea().maxY(), to, layer.flipped);
          drawDashedLine_(from, to, *painter);

          dataToWidget(distance_item->getEndPoint().getX(), 0, from, layer.flipped);

          dataToWidget(distance_item->getEndPoint().getX(), getVisibleArea().maxY(), to, layer.flipped);
          drawDashedLine_(from, to, *painter);
        }
      }
      QPoint begin, end; 
      switch (draw_modes_[i])
      {
      case DM_PEAKS:
        //---------------------DRAWING PEAKS---------------------

        for (SpectrumConstIteratorType it = vbegin; it != vend; ++it)
        {
          if (layer.filters.passes(spectrum, it - spectrum.begin()))
          {
            // use peak colors stored in the layer, if available
            if (layer.peak_colors_1d.size() == spectrum.size())
            {
              // find correct peak index
              Size peak_index = std::distance(spectrum.begin(), it);
              pen.setColor(layer.peak_colors_1d[peak_index]);
              painter->setPen(pen);
            }

            // Warn if non-empty peak color array present but size doesn't match number of peaks
            // This indicates a bug but we gracefuly just issue a warning
            if (!layer.peak_colors_1d.empty() &&
                layer.peak_colors_1d.size() < spectrum.size())
            {
              OPENMS_LOG_ERROR << "Peak color array size ("
                               << layer.peak_colors_1d.size()
                               << ") doesn't match number of peaks ("
                               << spectrum.size()
                               << ") in spectrum."
                               << endl;
            }
            dataToWidget(*it, end, layer.flipped);
            dataToWidget(it->getMZ(), 0.0f, begin, layer.flipped);
            // draw peak
            painter->drawLine(begin, end);
          }
        }
        break;

      case DM_CONNECTEDLINES:
      {
        //---------------------DRAWING CONNECTED LINES---------------------

        QPainterPath path;

        // connect peaks in visible area; (no clipping needed)
        bool first_point = true;
        for (SpectrumConstIteratorType it = vbegin; it != vend; it++)
        {
          dataToWidget(*it, begin, layer.flipped);

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
          dataToWidget(*(vbegin - 1), begin, layer.flipped);
          dataToWidget(*(vbegin), end, layer.flipped);
          painter->drawLine(begin, end);
        }

        // clipping on right side
        if (vend != spectrum.end() && vend != spectrum.begin())
        {
          dataToWidget(*(vend - 1), begin, layer.flipped);
          dataToWidget(*(vend), end, layer.flipped);
          painter->drawLine(begin, end);
        }
      }
      break;

      default:
        throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }

      // annotate interesting m/z's
      if (draw_interesting_MZs_) { drawMZAtInterestingPeaks_(i, *painter); }

      // draw all annotation items
      drawAnnotations(i, *painter);

      // draw a legend
      if (param_.getValue("show_legend").toBool())
      {
        const SpectrumType & spectrum = getLayer(i).getCurrentSpectrum();
        double xpos = getVisibleArea().maxX() - (getVisibleArea().maxX() - getVisibleArea().minX()) * 0.1;
        SpectrumConstIteratorType tmp  = max_element(spectrum.MZBegin(visible_area_.minX()), spectrum.MZEnd(xpos), PeakType::IntensityLess());
        if (tmp != spectrum.end())
        {
          PointType position(xpos, std::max<double>(tmp->getIntensity() - 100, tmp->getIntensity() * 0.8));
          Annotation1DPeakItem item = Annotation1DPeakItem(position, layer.getName().toQString(), QColor(layer.param.getValue("peak_color").toQString()));
          item.draw(this, *painter);
        }
      }
    }

    if (mirror_mode_)
    {
      painter->save();

      if (!show_alignment_)
      {
        // draw x-axis
        painter->setPen(Qt::black);
        painter->drawLine(0, height() / 2, width(), height() / 2);
      }
      else
      {
        drawAlignment(*painter);
        // two x-axes:
        painter->setPen(Qt::black);
        painter->drawLine(0, height() / 2 + 5, width(), height() / 2 + 5);
        painter->drawLine(0, height() / 2 - 5, width(), height() / 2 - 5);
      }
      painter->restore();
    }
    else   // !mirror_mode_
    {
      if (show_alignment_)
      {
        drawAlignment(*painter);
      }
    }

    // draw measuring line when in measure mode and valid measurement start peak selected
    if (action_mode_ == AM_MEASURE && measurement_start_.isValid())
    {
      QPoint measurement_end_point(last_mouse_pos_.x(), measurement_start_point_.y());
      // draw a complete temporary Annotation1DDistanceItem which includes the distance;
      // as an alternative to a simple line: painter->drawLine(measurement_start_point_, measurement_end_point);
      Annotation1DDistanceItem::PointType ps(widgetToData(measurement_start_point_, true));
      Annotation1DDistanceItem::PointType pe(widgetToData(measurement_end_point, true));
      Annotation1DDistanceItem(QString::number(pe.getX() - ps.getX(), 'f', 4), ps, pe).draw(this, *painter, false);
    }
    // draw highlighted measurement start peak and selected peak
    bool with_elongation = (action_mode_ == AM_MEASURE);
    drawHighlightedPeak_(getCurrentLayerIndex(), measurement_start_, *painter, with_elongation);
    drawHighlightedPeak_(getCurrentLayerIndex(), selected_peak_, *painter, with_elongation);

    //draw delta for measuring
    if (action_mode_ == AM_MEASURE && measurement_start_.isValid())
    {
      drawDeltas_(*painter, measurement_start_, selected_peak_);
    }
    else
    {
      drawCoordinates_(*painter, selected_peak_);
    }

    // draw text box (supporting HTML) on the right side of the canvas
    if (!text_box_content_.isEmpty())
    {
      painter->save();
      double w = text_box_content_.size().width();
      double h = text_box_content_.size().height();
      //draw text
      painter->setPen(Qt::black);
      painter->translate(width() - w - 2, 3);
      painter->fillRect(static_cast<int>(width() - w - 2),
                        3,
                        static_cast<int>(w),
                        static_cast<int>(h),
                        QColor(255, 255, 255, 200));
      text_box_content_.drawContents(painter);
      painter->restore();
    }

    if (show_timing_)
    {
      cout << "paint event took " << timer.elapsed() << " ms" << endl;
    }
  }

  void Spectrum1DCanvas::drawHighlightedPeak_(Size layer_index, const PeakIndex& peak, QPainter& painter, bool draw_elongation)
  {
    if (peak.isValid())
    {
      QPoint begin;
      const auto& spec = getLayer(layer_index).getCurrentSpectrum();
      if (peak.peak >= spec.size())
      {
        // somehow the peak is invalid. This happens from time to time and should be tracked down elsewhere
        // but its hard to reproduce (changing spectra in 1D view using arrow keys while hovering over the spectrum with the mouse?).
        return;
      }
      const ExperimentType::PeakType& sel = spec[peak.peak];

      painter.setPen(QPen(QColor(param_.getValue("highlighted_peak_color").toQString()), 2));

      updatePercentageFactor_(layer_index);

      dataToWidget(sel, begin, getLayer(layer_index).flipped);
      QPoint top_end(begin);

      bool layer_flipped = getLayer(layer_index).flipped;
      if (isMzToXAxis())
      {
        if (layer_flipped)
        {
          top_end.setY(height());
        }
        else
        {
          top_end.setY(0);
        }
      }
      else
      {
        if (!layer_flipped)
        {
          top_end.setX(width());
        }
        else // should not happen
        {
          top_end.setX(0);
        }
      }

      // paint the crosshair only for currently selected peaks of the current layer
      if (layer_index == getCurrentLayerIndex() && (peak == measurement_start_ || peak == selected_peak_))
      {
        painter.drawLine(begin.x(), begin.y() - 4, begin.x(), begin.y() + 4);
        painter.drawLine(begin.x() - 4, begin.y(), begin.x() + 4, begin.y());
      }
      // draw elongation as dashed line (while in measure mode and for all existing distance annotations)
      if (draw_elongation)
      {
        drawDashedLine_(begin, top_end, painter);
      }
    }
  }

  void Spectrum1DCanvas::drawDashedLine_(const QPoint& from, const QPoint& to, QPainter& painter)
  {
    QPen pen;
    QVector<qreal> dashes;
    dashes << 5 << 5 << 1 << 5;
    pen.setDashPattern(dashes);
    pen.setColor(QColor(param_.getValue("highlighted_peak_color").toQString()));
    painter.save();
    painter.setPen(pen);
    painter.drawLine(from, to);
    painter.restore();
  }

  void Spectrum1DCanvas::drawAnnotations(Size layer_index, QPainter& painter)
  {
    LayerData& layer = getLayer(layer_index);
    updatePercentageFactor_(layer_index);
    QColor col{ QColor(layer.param.getValue("annotation_color").toQString()) };
    // 0: default pen; 1: selected pen
    QPen pen[2] = { col, col.lighter() };

    // TODO: remove - just some debug code
    if (layer.getCurrentAnnotations().empty())
    {
      QColor col{ QColor(layer.param.getValue("annotation_color").toQString()) };
      for (double mz = 0.0; mz <= 2000.0; mz += 100.0)
      {
        auto item = new Annotation1DVerticalLineItem(mz, col, "Test");
        getCurrentLayer().getCurrentAnnotations().push_front(item);
      }
    }

    for (const auto& c : layer.getCurrentAnnotations())
    {
      painter.setPen(pen[c->isSelected()]);
      c->draw(this, painter, layer.flipped);
    }
  }

  void Spectrum1DCanvas::drawMZAtInterestingPeaks_(Size layer_index, QPainter& painter)
  {
    LayerData& layer = getLayer(layer_index);
    const MSSpectrum& current_spectrum = layer.getCurrentSpectrum();

    bool flipped = layer.flipped;
    updatePercentageFactor_(layer_index);

    // get visible peaks
    auto vbegin = current_spectrum.MZBegin(visible_area_.minX());
    auto vend = current_spectrum.MZEnd(visible_area_.maxX());

    if (vbegin == vend) { return; }

    // find interesting peaks

    // copy visible peaks into spec
    MSSpectrum spec;
    for (auto it(vbegin); it != vend; ++it) { spec.push_back(*it); }

    // calculate distance between first and last peak
    --vend;
    double visible_range = vend->getMZ() - vbegin->getMZ();

    // remove 0 intensities
    ThresholdMower threshold_mower_filter;
    threshold_mower_filter.filterPeakSpectrum(spec);

    // deisotope so we don't consider higher isotopic peaks
    Deisotoper::deisotopeAndSingleCharge(spec,
      100,     // tolerance
      true,   // ppm
      1, 6,   // min / max charge
      false,  // keep only deisotoped
      3, 10,  // min / max isopeaks
      false,  // don't convert fragment m/z to mono-charge
      true);  // annotate charge in integer data array

    // filter for local high-intensity peaks
    WindowMower window_mower_filter;
    Param filter_param = window_mower_filter.getParameters();
    double window_size = visible_range / 10.0;
    filter_param.setValue("windowsize", window_size, "The size of the sliding window along the m/z axis.");
    filter_param.setValue("peakcount", 2, "The number of peaks that should be kept.");
    filter_param.setValue("movetype", "slide", "Whether sliding window (one peak steps) or jumping window (window size steps) should be used.");
    window_mower_filter.setParameters(filter_param);
    window_mower_filter.filterPeakSpectrum(spec);

    NLargest nlargest_filter(10);  // maximum number of annotated m/z's in visible area
    nlargest_filter.filterPeakSpectrum(spec);
    spec.sortByPosition(); // nlargest changes order

    for (Size i = 0; i != spec.size(); ++i)
    {
      Size current_peak_index = current_spectrum.findNearest(spec[i].getMZ());
      double mz(current_spectrum[current_peak_index].getMZ());
      double intensity(current_spectrum[current_peak_index].getIntensity());

      QString label = String::number(mz, 4).toQString();

      if (!spec.getIntegerDataArrays().empty()
          && spec.getIntegerDataArrays()[0].size() == spec.size())
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
      item.draw(this, painter, flipped);
    }
  }

  void Spectrum1DCanvas::changeVisibleArea_(const AreaType& new_area, bool repaint, bool add_to_stack)
  {
    // set new visible area (if changed)
    if (new_area != visible_area_)
    {
      visible_area_ = new_area;
      updateScrollbars_();
      recalculateSnapFactor_();
      emit visibleAreaChanged(new_area);
    }

    // store old zoom state
    if (add_to_stack) { zoomAdd_(new_area); }

    // repaint
    if (repaint) { update_(OPENMS_PRETTY_FUNCTION); }
  }

  bool Spectrum1DCanvas::finishAdding_()
  {
    if (getCurrentLayer().type != LayerData::DT_PEAK)
    {
      QMessageBox::critical(this, "Error", "This widget supports peak data only. Aborting!");
      return false;
    }

    getCurrentLayer().updateRanges();

    // Abort if no data points are contained (note that all data could be on disk)
    if (getCurrentLayer().getCurrentSpectrum().empty())
    {
      popIncompleteLayer_("Cannot add a dataset that contains no survey scans. Aborting!");
      return false;
    }

    const MSSpectrum& spectrum = getCurrentLayer().getCurrentSpectrum();

    // add new draw mode and style (default: peaks)
    draw_modes_.push_back(DM_PEAKS);
    SpectrumSettings::SpectrumType spectrum_type = spectrum.getType(true);

    if (spectrum_type == SpectrumSettings::PROFILE)
    {
      draw_modes_.back() = DM_CONNECTEDLINES;
    }
    peak_penstyle_.push_back(Qt::SolidLine);


    // Change peak color if this is not the first layer
    switch (getCurrentLayerIndex() % 5)
    {
    case 0:
      getCurrentLayer().param.setValue("peak_color", "#0000ff");
      getCurrentLayer().param.setValue("annotation_color", "#005500");
      break;

    case 1:
      getCurrentLayer().param.setValue("peak_color", "#00cc00");
      getCurrentLayer().param.setValue("annotation_color", "#005500");
      break;

    case 2:
      getCurrentLayer().param.setValue("peak_color", "#cc0000");
      getCurrentLayer().param.setValue("annotation_color", "#550055");
      break;

    case 3:
      getCurrentLayer().param.setValue("peak_color", "#00cccc");
      getCurrentLayer().param.setValue("annotation_color", "#005555");
      break;

    case 4:
      getCurrentLayer().param.setValue("peak_color", "#ffaa00");
      getCurrentLayer().param.setValue("annotation_color", "#550000");
      break;
    }

    // sort spectra in ascending order of position (ensure that we sort all spectra as well as the currently
    // TODO: check why this is need since we load data already sorted! 
    for (Size i = 0; i < getCurrentLayer().getPeakData()->size(); ++i)
    {
      (*getCurrentLayer().getPeakDataMuteable())[i].sortByPosition();
    }
    getCurrentLayer().sortCurrentSpectrumByPosition();

    getCurrentLayer().annotations_1d.resize(getCurrentLayer().getPeakData()->size());

    // update nearest peak
    selected_peak_.clear();

    // update ranges
    recalculateRanges_(0, 2, 1);
    double width = overall_data_range_.width();
    overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
    overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
    overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());
    resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway

    // warn if negative intensities are contained
    if (getCurrentMinIntensity() < 0.0)
    {
      QMessageBox::warning(this, "Warning", "This dataset contains negative intensities. Use it at your own risk!");
    }

    if (getLayerCount() == 2)
    {
      setIntensityMode(IM_PERCENTAGE);
    }

    emit layerActivated(this);

    return true;
  }

  void Spectrum1DCanvas::drawCoordinates_(QPainter& painter, const PeakIndex& peak)
  {
    if (!peak.isValid())
      return;

    //determine coordinates;
    double mz = 0.0;
    float it = 0.0;
    // only peak data is supported here
    if (getCurrentLayer().type != LayerData::DT_PEAK)
    {
      QMessageBox::critical(this, "Error", "This widget supports peak data only. Aborting!");
      return;
    }
    mz = getCurrentLayer().getCurrentSpectrum()[peak.peak].getMZ();
    it = getCurrentLayer().getCurrentSpectrum()[peak.peak].getIntensity();

    //draw text
    QStringList lines;
    String text;
    int precision(2);

    if (isMzToXAxis() ^ is_swapped_) // XOR
    { // only if either one of the conditions holds
      text = "RT:  "; // two spaces, ensuring same indentation as "m/z: " and "int: "
      precision = 2;
    }
    else // only if none or both are true
    {
      text = "m/z: ";
      precision = 8;
    }
    lines.push_back(text.c_str() +  QLocale::c().toString(mz, 'f', precision));  // adds group separators (consistency with intensity)
    lines.push_back("Int: " + QLocale::c().toString(it, 'f', 2));                // adds group separators (every 1e3), to better visualize large numbers (e.g. 23.009.646.54,3));
    drawText_(painter, lines);
  }

  void Spectrum1DCanvas::drawDeltas_(QPainter& painter, const PeakIndex& start, const PeakIndex& end)
  {
    if (!start.isValid())
      return;

    //determine coordinates;
    double mz;
    float it;
    float ppm;

    if (getCurrentLayer().type != LayerData::DT_PEAK)
    {
      QMessageBox::critical(this, "Error", "This widget supports peak data only. Aborting!");
      return;
    }

    if (end.isValid())
    {
      mz = getCurrentLayer().getCurrentSpectrum()[end.peak].getMZ() - getCurrentLayer().getCurrentSpectrum()[start.peak].getMZ();
      it = getCurrentLayer().getCurrentSpectrum()[end.peak].getIntensity() - getCurrentLayer().getCurrentSpectrum()[start.peak].getIntensity();
    }
    else
    {
      PointType point = widgetToData_(last_mouse_pos_);
      mz = point[0] - getCurrentLayer().getCurrentSpectrum()[start.peak].getMZ();
      it = std::numeric_limits<double>::quiet_NaN();
    }
    ppm = (mz / getCurrentLayer().getCurrentSpectrum()[start.peak].getMZ()) * 1e6;

    //draw text
    QStringList lines;
    String text;
    int precision(2);
    if (isMzToXAxis() ^ is_swapped_) // XOR
    { // only if either one of the conditions holds
      text = "RT delta: ";
      precision = 2;
    }
    else // only if none or both are true
    {
      text = "m/z delta: ";
      precision = 6;
    }
    lines.push_back(text.c_str() + QString::number(mz, 'f', precision) + " (" + QString::number(ppm, 'f', 1) +" ppm)");

    if (boost::math::isinf(it) || boost::math::isnan(it))
    {
      lines.push_back("Int ratio: n/a");
    }
    else
    {
      lines.push_back("Int ratio: " + QString::number(it, 'f', 2));
    }
    drawText_(painter, lines);
  }

  void Spectrum1DCanvas::recalculateSnapFactor_()
  {
    if (intensity_mode_ == IM_SNAP)
    {
      double local_max  = -numeric_limits<double>::max();
      for (Size i = 0; i < getLayerCount(); ++i)
      {
        const SpectrumType & spectrum = getLayer(i).getCurrentSpectrum();
        SpectrumConstIteratorType tmp = max_element(spectrum.MZBegin(visible_area_.minX()), spectrum.MZEnd(visible_area_.maxX()), PeakType::IntensityLess());
        if (tmp != spectrum.end() && tmp->getIntensity() > local_max)
        {
          local_max = tmp->getIntensity();
        }
      }

      // add some margin on top of local maximum to be sure we are able to draw labels inside the view
      snap_factors_[0] = overall_data_range_.maxPosition()[1] / (local_max * TOP_MARGIN);
    }
    else if (intensity_mode_ == IM_PERCENTAGE)
    {
      snap_factors_[0] = 1.0 / TOP_MARGIN;
    }
    else
    {
      snap_factors_[0] = 1.0;
    }
  }

  void Spectrum1DCanvas::updateScrollbars_()
  {
    emit updateHScrollbar(overall_data_range_.minPosition()[0], visible_area_.minPosition()[0], visible_area_.maxPosition()[0], overall_data_range_.maxPosition()[0]);
    emit updateVScrollbar(1, 1, 1, 1);
  }

  void Spectrum1DCanvas::horizontalScrollBarChange(int value)
  {
    changeVisibleArea_(value, value + (visible_area_.maxPosition()[0] - visible_area_.minPosition()[0]));
  }

  void Spectrum1DCanvas::showCurrentLayerPreferences()
  {
    Internal::Spectrum1DPrefDialog dlg(this);
    LayerData& layer = getCurrentLayer();

    ColorSelector* peak_color = dlg.findChild<ColorSelector*>("peak_color");
    ColorSelector* icon_color = dlg.findChild<ColorSelector*>("icon_color");
    ColorSelector* annotation_color = dlg.findChild<ColorSelector*>("annotation_color");
    ColorSelector* bg_color = dlg.findChild<ColorSelector*>("bg_color");
    ColorSelector* selected_color = dlg.findChild<ColorSelector*>("selected_color");

    peak_color->setColor(QColor(layer.param.getValue("peak_color").toQString()));
    icon_color->setColor(QColor(layer.param.getValue("icon_color").toQString()));
    annotation_color->setColor(QColor(layer.param.getValue("annotation_color").toQString()));
    bg_color->setColor(QColor(param_.getValue("background_color").toQString()));
    selected_color->setColor(QColor(param_.getValue("highlighted_peak_color").toQString()));

    if (dlg.exec())
    {
      layer.param.setValue("peak_color", peak_color->getColor().name());
      layer.param.setValue("icon_color", icon_color->getColor().name());
      layer.param.setValue("annotation_color", annotation_color->getColor().name());
      param_.setValue("background_color", bg_color->getColor().name());
      param_.setValue("highlighted_peak_color", selected_color->getColor().name());

      emit preferencesChange();
    }
  }

  void Spectrum1DCanvas::currentLayerParamtersChanged_()
  {
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum1DCanvas::contextMenuEvent(QContextMenuEvent* e)
  {
    if (layers_.empty()) { return; }

    QMenu* context_menu = new QMenu(this);
    QAction* result = nullptr;

    Annotations1DContainer& annots_1d = getCurrentLayer().getCurrentAnnotations();
    Annotation1DItem* annot_item = annots_1d.getItemAt(e->pos());
    if (annot_item)
    {
      annots_1d.deselectAll();
      annots_1d.selectItemAt(e->pos());
      update_(OPENMS_PRETTY_FUNCTION);

      context_menu->addAction("Edit");
      context_menu->addAction("Delete");
      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->text() == "Delete")
        {
          vector<Annotation1DItem*> as;
          as.push_back(annot_item);
          getCurrentLayer().removePeakAnnotationsFromPeptideHit(as);
          annots_1d.removeSelectedItems();
        }
        else if (result->text() == "Edit")
        {
          annot_item->editText();
          getCurrentLayer().synchronizePeakAnnotations();
        }
        update_(OPENMS_PRETTY_FUNCTION);
      }
    }
    else // !annot_item
    {
      //Display name and warn if current layer invisible
      String layer_name = String("Layer: ") + getCurrentLayer().getName();
      if (!getCurrentLayer().visible)
      {
        layer_name += " (invisible)";
      }
      context_menu->addAction(layer_name.toQString())->setEnabled(false);
      context_menu->addSeparator();
      QAction* new_action = context_menu->addAction("Add label");
      if (mirror_mode_ && (getCurrentLayer().flipped ^ (e->pos().y() > height() / 2)))
      {
        new_action->setEnabled(false);
      }
      new_action = context_menu->addAction("Add peak annotation");
      PeakIndex near_peak = findPeakAtPosition_(e->pos());
      if (!near_peak.isValid())
      {
        new_action->setEnabled(false);
      }
      new_action = context_menu->addAction("Add peak annotation mz");
      if (!near_peak.isValid())
      {
        new_action->setEnabled(false);
      }
      context_menu->addSeparator();
      new_action = context_menu->addAction("Reset alignment");
      if (!show_alignment_)
      {
        new_action->setEnabled(false);
      }
      context_menu->addSeparator();

      context_menu->addAction("Layer meta data");

      QMenu* save_menu = new QMenu("Save");
      save_menu->addAction("Layer");
      save_menu->addAction("Visible layer data");
      save_menu->addAction("As image");

      QMenu* settings_menu = new QMenu("Settings");
      settings_menu->addAction("Show/hide grid lines");
      settings_menu->addAction("Show/hide axis legends");
      settings_menu->addAction("Style: Stick <--> Area");
      settings_menu->addAction("Intensity: Absolute <--> Percent");
      settings_menu->addAction("Show/hide ion ladder in ID view");
      settings_menu->addAction("Show/hide automated m/z annotations");
      settings_menu->addSeparator();
      settings_menu->addAction("Preferences");

      context_menu->addMenu(save_menu);
      context_menu->addMenu(settings_menu);

      // only add to context menu if there is a MS1 map
      if (getCurrentLayer().getPeakData()->containsScanOfLevel(1))
      {
        context_menu->addAction("Switch to 2D view");
        context_menu->addAction("Switch to 3D view");
      }

      if (getCurrentLayer().getCurrentSpectrum().containsIMData())
      {
        context_menu->addAction("Switch to ion mobility view");
      }

      if (getCurrentLayer().isDIAData())
      {
        context_menu->addAction("Switch to DIA-MS view");
      }

      //add external context menu
      if (context_add_)
      {
        context_menu->addSeparator();
        context_menu->addMenu(context_add_);
      }

      //evaluate menu
      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->text() == "Preferences")
        {
          showCurrentLayerPreferences();
        }
        else if (result->text() == "Show/hide grid lines")
        {
          showGridLines(!gridLinesShown());
        }
        else if (result->text() == "Show/hide axis legends")
        {
          emit changeLegendVisibility();
        }
        else if (result->text() == "Show/hide automated m/z annotations")
        {
          setDrawInterestingMZs(!draw_interesting_MZs_);
        }
        else if (result->text() == "Layer" || result->text() == "Visible layer data")
        {
          saveCurrentLayer(result->text() == "Visible layer data");
        }
        else if (result->text() == "As image")
        {
          spectrum_widget_->saveAsImage();
        }
        else if (result->text() == "Style: Stick <--> Area")
        {
          if (getDrawMode() != DM_PEAKS)
          {
            setDrawMode(DM_PEAKS);
          }
          else
          {
            setDrawMode(DM_CONNECTEDLINES);
          }
        }
        else if (result->text() == "Intensity: Absolute <--> Percent")
        {
          if (getIntensityMode() != IM_PERCENTAGE)
          {
            setIntensityMode(IM_PERCENTAGE);
          }
          else
          {
            setIntensityMode(IM_SNAP);
          }
        }
        else if (result->text() == "Layer meta data")
        {
          showMetaData(true);
        }
        else if (result->text() == "Add label")
        {
          addUserLabelAnnotation_(e->pos());
        }
        else if (result->text() == "Add peak annotation")
        {
          addUserPeakAnnotation_(near_peak);
        }
        else if (result->text() == "Add peak annotation mz")
        {
          QString label = String::number(getCurrentLayer().getCurrentSpectrum()[near_peak.peak].getMZ(), 4).toQString();
          addPeakAnnotation(near_peak, label, getCurrentLayer().param.getValue("peak_color").toQString());
        }
        else if (result->text() == "Reset alignment")
        {
          resetAlignment();
        }
        else if (result->text() == "Switch to 2D view")
        {
          emit showCurrentPeaksAs2D();
        }
        else if (result->text() == "Switch to 3D view")
        {
          emit showCurrentPeaksAs3D();
        }
        else if (result->text() == "Switch to ion mobility view")
        {
          emit showCurrentPeaksAsIonMobility();
        }
        else if (result->text() == "Switch to DIA-MS view")
        {
          emit showCurrentPeaksAsDIA();
        }
        else if (result->text() == "Show/hide ion ladder in ID view")
        {
          // toggle visibility of ion ladder
          setIonLadderVisible(!isIonLadderVisible());
        }
      }
    }
    e->accept();
  }

  void Spectrum1DCanvas::setTextBox(const QString& html)
  {
    text_box_content_.setHtml(html);
  }

  void Spectrum1DCanvas::addUserLabelAnnotation_(const QPoint& screen_position)
  {
    bool ok;
    QString text = QInputDialog::getText(this, "Add label", "Enter text:", QLineEdit::Normal, "", &ok);
    if (ok && !text.isEmpty())
    {
      addLabelAnnotation_(screen_position, text);
    }
  }

  void Spectrum1DCanvas::addLabelAnnotation_(const QPoint& screen_position, QString text)
  {
    updatePercentageFactor_(getCurrentLayerIndex());

    PointType position = widgetToData(screen_position, true);
    Annotation1DItem* item = new Annotation1DTextItem(position, text);
    getCurrentLayer().getCurrentAnnotations().push_front(item);

    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum1DCanvas::addUserPeakAnnotation_(PeakIndex near_peak)
  {
    bool ok;
    QString text = QInputDialog::getText(this, "Add peak annotation", "Enter text:", QLineEdit::Normal, "", &ok);
    if (ok && !text.isEmpty())
    {
      addPeakAnnotation(near_peak, text, QColor(getCurrentLayer().param.getValue("peak_color").toQString()));
    }
  }

  Annotation1DItem* Spectrum1DCanvas::addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color)
  {
    PeakType peak = getCurrentLayer().getCurrentSpectrum()[peak_index.peak];
    PointType position(peak.getMZ(), peak.getIntensity());
    Annotation1DItem* item = new Annotation1DPeakItem(position, text, color);
    item->setSelected(false);
    getCurrentLayer().getCurrentAnnotations().push_front(item);
    update_(OPENMS_PRETTY_FUNCTION);
    return item;
  }

  void Spectrum1DCanvas::saveCurrentLayer(bool visible)
  {
    const LayerData& layer = getCurrentLayer();

    //determine proposed filename
    String proposed_name = param_.getValue("default_path");
    if (!visible && !layer.filename.empty())
    {
      proposed_name = layer.filename;
    }

    QString selected_filter = "";
    QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(), "mzML files (*.mzML);;mzData files (*.mzData);;mzXML files (*.mzXML);;All files (*)", &selected_filter);
    if (!file_name.isEmpty())
    {
      // check whether a file type suffix has been given
      // first check mzData and mzXML then mzML
      // if the setting is at "All files"
      // mzML will be used
      String upper_filename = file_name;
      upper_filename.toUpper();
      if (selected_filter == "mzData files (*.mzData)")
      {
        if (!upper_filename.hasSuffix(".MZDATA"))
        {
          file_name += ".mzData";
        }
      }
      else if (selected_filter == "mzXML files (*.mzXML)")
      {
        if (!upper_filename.hasSuffix(".MZXML"))
        {
          file_name += ".mzXML";
        }
      }
      else
      {
        if (!upper_filename.hasSuffix(".MZML"))
        {
          file_name += ".mzML";
        }
      }

      if (visible)
      {
        ExperimentType out;
        getVisiblePeakData(out);
        addDataProcessing_(out, DataProcessing::FILTERING);
        FileHandler().storeExperiment(file_name, out, ProgressLogger::GUI);
      }
      else
      {
        // TODO: this will not work if the data is cached on disk
        FileHandler().storeExperiment(file_name, *layer.getPeakData(), ProgressLogger::GUI);
      }
    }
  }

  bool Spectrum1DCanvas::flippedLayersExist()
  {
    for (Size i = 0; i < getLayerCount(); ++i)
    {
      if (layers_.getLayer(i).flipped)
      {
        return true;
      }
    }
    return false;
  }

  void Spectrum1DCanvas::updateLayer(Size i)
  {
    //update nearest peak
    selected_peak_.clear();

    //update ranges
    recalculateRanges_(0, 2, 1);
    double width = overall_data_range_.width();
    overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
    overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
    overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());

    resetZoom();
    modificationStatus_(i, false);
  }

  void Spectrum1DCanvas::zoom_(int x, int y, bool zoom_in)
  {
    if (!zoom_in)
    {
      zoomBack_();
    }
    else
    {
      const PointType::CoordinateType zoom_factor = 0.8;
      AreaType new_area;
      if (isMzToXAxis())
      {
        new_area.setMinX(visible_area_.min_[0] + (1.0 - zoom_factor) * (visible_area_.max_[0] - visible_area_.min_[0]) * (PointType::CoordinateType)x / width());
        new_area.setMaxX(new_area.min_[0] + zoom_factor * (visible_area_.max_[0] - visible_area_.min_[0]));
        new_area.setMinY(visible_area_.minY());
        new_area.setMaxY(visible_area_.maxY());
      }
      else
      {
        new_area.setMinX(visible_area_.min_[0] + (1.0 - zoom_factor) * (visible_area_.max_[0] - visible_area_.min_[0]) * (PointType::CoordinateType)(height() - y) / height());
        new_area.setMaxX(new_area.min_[0] + zoom_factor * (visible_area_.max_[0] - visible_area_.min_[0]));
        new_area.setMinY(visible_area_.minY());
        new_area.setMaxY(visible_area_.maxY());
      }

      if (new_area != visible_area_)
      {
        zoomAdd_(new_area);
        zoom_pos_ = --zoom_stack_.end(); // set to last position
        changeVisibleArea_(*zoom_pos_);
      }
    }
  }

  /// Go forward in zoom history
  void Spectrum1DCanvas::zoomForward_()
  {
    // if at end of zoom level then simply add a new zoom
    if (zoom_pos_ == zoom_stack_.end() || (zoom_pos_ + 1) == zoom_stack_.end())
    {
      AreaType new_area;
      // distance of areas center to border times a zoom factor of 0.8
      AreaType::CoordinateType size0 = visible_area_.width() / 2 * 0.8;
      new_area.setMinX(visible_area_.center()[0] - size0);
      new_area.setMaxX(visible_area_.center()[0] + size0);
      new_area.setMinY(visible_area_.minY());
      new_area.setMaxY(visible_area_.maxY());
      zoomAdd_(new_area);
      zoom_pos_ = --zoom_stack_.end(); // set to last position
    }
    else
    {
      // goto next zoom level
      ++zoom_pos_;
    }
    changeVisibleArea_(*zoom_pos_);
  }

  void Spectrum1DCanvas::translateLeft_(Qt::KeyboardModifiers m)
  {
    double newLo = visible_area_.minX();
    double newHi = visible_area_.maxX();
    if (m == Qt::NoModifier)
    { // 5% shift
      double shift = 0.05 * visible_area_.width();
      newLo -= shift;
      newHi -= shift;
    }
    else if (m == Qt::ShiftModifier)
    { // jump to the next peak (useful for sparse data)
      const LayerData::ExperimentType::SpectrumType& spec = getCurrentLayer().getCurrentSpectrum();
      PeakType p_temp(visible_area_.minX(), 0);
      SpectrumConstIteratorType it_next = lower_bound(spec.begin(), spec.end(), p_temp, PeakType::MZLess()); // find first peak in current range
      if (it_next != spec.begin()) --it_next; // move one peak left
      if (it_next == spec.end()) return;
      newLo = it_next->getMZ() - visible_area_.width() / 2; // center the next peak to the left
      newHi = it_next->getMZ() + visible_area_.width() / 2;
    }

    // check if we are falling out of bounds
    if (newLo < overall_data_range_.minX())
    {
      newLo = overall_data_range_.minX();
      newHi = newLo + visible_area_.width();
    }
    // change data area
    changeVisibleArea_(newLo, newHi);
  }

  void Spectrum1DCanvas::translateRight_(Qt::KeyboardModifiers m)
  {
    double newLo = visible_area_.minX();
    double newHi = visible_area_.maxX();
    if (m == Qt::NoModifier)
    { // 5% shift
      double shift = 0.05 * visible_area_.width();
      newLo += shift;
      newHi += shift;
    }
    else if (m == Qt::ShiftModifier)
    { // jump to the next peak (useful for sparse data)
      const LayerData::ExperimentType::SpectrumType& spec = getCurrentLayer().getCurrentSpectrum();
      PeakType p_temp(visible_area_.maxX(), 0);
      SpectrumConstIteratorType it_next = upper_bound(spec.begin(), spec.end(), p_temp, PeakType::MZLess()); // first right-sided peak outside the current range
      if (it_next == spec.end()) return;
      newLo = it_next->getMZ() - visible_area_.width() / 2; // center the next peak to the right
      newHi = it_next->getMZ() + visible_area_.width() / 2;
    }

    // check if we are falling out of bounds
    if (newHi > overall_data_range_.maxX())
    {
      newHi = overall_data_range_.maxX();
      newLo = newHi - visible_area_.width();
    }
    // change data area
    changeVisibleArea_(newLo, newHi);
  }

  /// Returns whether this widget is currently in mirror mode
  bool Spectrum1DCanvas::mirrorModeActive()
  {
    return mirror_mode_;
  }

  /// Sets whether this widget is currently in mirror mode
  void Spectrum1DCanvas::setMirrorModeActive(bool b)
  {
    mirror_mode_ = b;
    qobject_cast<Spectrum1DWidget*>(spectrum_widget_)->toggleMirrorView(b);
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum1DCanvas::paintGridLines_(QPainter& painter)
  {
    if (!show_grid_ || !spectrum_widget_)
      return;

    QPen p1(QColor(130, 130, 130));
    p1.setStyle(Qt::DashLine);
    QPen p2(QColor(170, 170, 170));
    p2.setStyle(Qt::DotLine);

    painter.save();

    unsigned int xl, xh, yl, yh;     //width/height of the diagram area, x, y coordinates of lo/hi x,y values

    xl = 0;
    xh = width();

    yl = height();
    yh = 0;

    // drawing of grid lines and associated text
    for (Size j = 0; j != spectrum_widget_->xAxis()->gridLines().size(); j++)
    {
      // style definitions
      switch (j)
      {
      case 0:           // style settings for big intervals
        painter.setPen(p1);
        break;

      case 1:           // style settings for small intervals
        painter.setPen(p2);
        break;

      default:
        std::cout << "empty vertical grid line vector error!" << std::endl;
        painter.setPen(QPen(QColor(0, 0, 0)));
        break;
      }

      for (std::vector<double>::const_iterator it = spectrum_widget_->xAxis()->gridLines()[j].begin(); it != spectrum_widget_->xAxis()->gridLines()[j].end(); ++it)
      {
        int x = static_cast<int>(Math::intervalTransformation(*it, spectrum_widget_->xAxis()->getAxisMinimum(), spectrum_widget_->xAxis()->getAxisMaximum(), xl, xh));
        painter.drawLine(x, yl, x, yh);
      }
    }

    for (Size j = 0; j != spectrum_widget_->yAxis()->gridLines().size(); j++)
    {

      // style definitions
      switch (j)
      {
      case 0:           // style settings for big intervals
        painter.setPen(p1);
        break;

      case 1:           // style settings for small intervals
        painter.setPen(p2);
        break;

      default:
        std::cout << "empty vertical grid line vector error!" << std::endl;
        painter.setPen(QPen(QColor(0, 0, 0)));
        break;
      }

      for (std::vector<double>::const_iterator it = spectrum_widget_->yAxis()->gridLines()[j].begin(); it != spectrum_widget_->yAxis()->gridLines()[j].end(); ++it)
      {
        int y = static_cast<int>(Math::intervalTransformation(*it, spectrum_widget_->yAxis()->getAxisMinimum(), spectrum_widget_->yAxis()->getAxisMaximum(), yl, yh));
        if (!mirror_mode_)
        {
          painter.drawLine(xl, y, xh, y);
        }
        else
        {
          if (!show_alignment_)
          {
            painter.drawLine(xl, y / 2, xh, y / 2);
            painter.drawLine(xl, yl - y / 2, xh, yl - y / 2);
          }
          else
          {
            double alignment_shrink_factor = 1.0;
            if (height() > 10)
            {
              alignment_shrink_factor = (double)(height() - 10) / (double)height();
            }
            painter.drawLine(xl, (int)((double)(y) * alignment_shrink_factor / 2.0), xh, (int)((double)(y) * alignment_shrink_factor / 2.0));
            painter.drawLine(xl, yl - (int)((double)(y) * alignment_shrink_factor / 2.0), xh, yl - (int)((double)(y) * alignment_shrink_factor / 2.0));
          }
        }
      }
    }

    painter.restore();
  }

  void Spectrum1DCanvas::performAlignment(Size layer_index_1, Size layer_index_2, const Param& param)
  {
    alignment_layer_1_ = layer_index_1;
    alignment_layer_2_ = layer_index_2;
    aligned_peaks_mz_delta_.clear();
    aligned_peaks_indices_.clear();

    if (layer_index_1 >= getLayerCount() || layer_index_2 >= getLayerCount())
    {
      return;
    }
    LayerData& layer_1 = getLayer(layer_index_1);
    LayerData& layer_2 = getLayer(layer_index_2);
    const ExperimentType::SpectrumType& spectrum_1 = layer_1.getCurrentSpectrum();
    const ExperimentType::SpectrumType& spectrum_2 = layer_2.getCurrentSpectrum();

    SpectrumAlignment aligner;
    aligner.setParameters(param);
    aligner.getSpectrumAlignment(aligned_peaks_indices_, spectrum_1, spectrum_2);

    for (Size i = 0; i < aligned_peaks_indices_.size(); ++i)
    {
      double line_begin_mz = spectrum_1[aligned_peaks_indices_[i].first].getMZ();
      double line_end_mz = spectrum_2[aligned_peaks_indices_[i].second].getMZ();
      aligned_peaks_mz_delta_.emplace_back(line_begin_mz, line_end_mz);
    }

    show_alignment_ = true;
    update_(OPENMS_PRETTY_FUNCTION);

    SpectrumAlignmentScore scorer;
    scorer.setParameters(param);

    alignment_score_ = scorer(spectrum_1, spectrum_2);
  }

  void Spectrum1DCanvas::resetAlignment()
  {
    aligned_peaks_indices_.clear();
    aligned_peaks_mz_delta_.clear();
    qobject_cast<Spectrum1DWidget*>(spectrum_widget_)->resetAlignment();
    show_alignment_ = false;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum1DCanvas::drawAlignment(QPainter& painter)
  {
    painter.save();

    //draw peak-connecting lines between the two spectra
    painter.setPen(Qt::red);
    QPoint begin_p, end_p;
    if (mirror_mode_)
    {
      double dummy = 0.0;
      for (Size i = 0; i < getAlignmentSize(); ++i)
      {
        dataToWidget(aligned_peaks_mz_delta_[i].first, dummy, begin_p);
        dataToWidget(aligned_peaks_mz_delta_[i].second, dummy, end_p);
        painter.drawLine(begin_p.x(), height() / 2 - 5, end_p.x(), height() / 2 + 5);
      }
    }
    else
    {
      const ExperimentType::SpectrumType& spectrum_1 = getLayer(alignment_layer_1_).getCurrentSpectrum();
      updatePercentageFactor_(alignment_layer_1_);
      for (Size i = 0; i < getAlignmentSize(); ++i)
      {
        dataToWidget(spectrum_1[aligned_peaks_indices_[i].first].getMZ(), 0, begin_p, false, true);
        dataToWidget(spectrum_1[aligned_peaks_indices_[i].first].getMZ(), spectrum_1[aligned_peaks_indices_[i].first].getIntensity(), end_p, false, true);
        painter.drawLine(begin_p.x(), begin_p.y(), end_p.x(), end_p.y());
      }
    }
    painter.restore();
  }

  Size Spectrum1DCanvas::getAlignmentSize()
  {
    return aligned_peaks_mz_delta_.size();
  }

  double Spectrum1DCanvas::getAlignmentScore()
  {
    return alignment_score_;
  }

  void Spectrum1DCanvas::intensityModeChange_()
  {
    recalculateSnapFactor_();
    ensureAnnotationsWithinDataRange_();
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Spectrum1DCanvas::ensureAnnotationsWithinDataRange_()
  {
    for (Size i = 0; i < getLayerCount(); ++i)
    {
      updatePercentageFactor_(i);
      Annotations1DContainer& ann_1d = getLayer(i).getCurrentAnnotations();
      for (Annotations1DContainer::Iterator it = ann_1d.begin(); it != ann_1d.end(); ++it)
      {
        (*it)->ensureWithinDataRange(this);
      }
    }
  }

  void Spectrum1DCanvas::updatePercentageFactor_(Size layer_index)
  {
    if (intensity_mode_ == IM_PERCENTAGE)
    {
      percentage_factor_ = overall_data_range_.maxPosition()[1] / getLayer(layer_index).getCurrentSpectrum().getMaxInt();
    }
    else
    {
      percentage_factor_ = 1.0;
    }
  }

  void Spectrum1DCanvas::flipLayer(Size index)
  {
    if (index < getLayerCount())
    {
      getLayer(index).flipped = !getLayer(index).flipped;
    }
  }

  void Spectrum1DCanvas::activateSpectrum(Size index, bool repaint)
  {
    // Note: even though the current spectrum may be on disk, there will still
    // be an in-memory representation in the peak data structure. Using
    // setCurrentSpectrumIndex will select the appropriate spectrum and load it
    // into memory.
    if (index < getCurrentLayer().getPeakData()->size())
    {
      getCurrentLayer().setCurrentSpectrumIndex(index);
      recalculateSnapFactor_();
      if (repaint)
      {
        update_(OPENMS_PRETTY_FUNCTION);
      }
    }
  }

  void Spectrum1DCanvas::setSwappedAxis(bool swapped)
  {
    is_swapped_ = swapped;
  }

  void Spectrum1DCanvas::setCurrentLayerPeakPenStyle(Qt::PenStyle ps)
  {
    // no layers
    if (layers_.empty())
    {
      return;
    }

    if (peak_penstyle_[getCurrentLayerIndex()] != ps)
    {
      peak_penstyle_[getCurrentLayerIndex()] = ps;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  std::vector<std::pair<Size, Size> > Spectrum1DCanvas::getAlignedPeaksIndices()
  {
    return aligned_peaks_indices_;
  }

  void Spectrum1DCanvas::setIonLadderVisible(bool show)
  {
    if (ion_ladder_visible_ != show)
    {
      ion_ladder_visible_ = show;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  void Spectrum1DCanvas::setDrawInterestingMZs(bool enable)
  {
    if (draw_interesting_MZs_ != enable)
    {
      draw_interesting_MZs_ = enable;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  bool Spectrum1DCanvas::isIonLadderVisible() const
  {
    return ion_ladder_visible_;
  }

  bool Spectrum1DCanvas::isDrawInterestingMZs() const
  {
    return draw_interesting_MZs_;
  }

} //Namespace
