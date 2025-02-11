// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/Plot1DCanvas.h>

// OpenMS
#include <OpenMS/VISUAL/PlotWidget.h>
#include <OpenMS/VISUAL/Plot1DWidget.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/DIALOGS/Plot1DPrefDialog.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotations1DContainer.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/RAIICleanup.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/COMPARISON/SpectrumAlignmentScore.h>
#include <OpenMS/COMPARISON/SpectrumAlignment.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <OpenMS/VISUAL/LayerData1DPeak.h>
#include <OpenMS/VISUAL/LayerData1DChrom.h>

// Qt
#include <QElapsedTimer>
#include <QMouseEvent>
#include <QPainter>
#include <QtWidgets/QInputDialog>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMessageBox>
#include <utility>


using namespace std;

namespace OpenMS
{
  using namespace Math;
  using namespace Internal;

  /// returns an MSExp with a single spec (converted from @p exp_sptr's chromatograms at index  @p index (or ondisc_sptr, if that should be empty)
  Plot1DCanvas::ExperimentSharedPtrType prepareChromatogram(Size index, const Plot1DCanvas::ExperimentSharedPtrType& exp_sptr, const Plot1DCanvas::ODExperimentSharedPtrType& ondisc_sptr)
  {
    // create a managed pointer fill it with a spectrum containing the chromatographic data
    LayerDataBase::ExperimentSharedPtrType chrom_exp_sptr(new LayerDataBase::ExperimentType());
    chrom_exp_sptr->setMetaValue("is_chromatogram", "true"); //this is a hack to store that we have chromatogram data
    LayerDataBase::ExperimentType::SpectrumType spectrum;

    // retrieve chromatogram (either from in-memory or on-disc representation)
    MSChromatogram current_chrom = exp_sptr->getChromatograms()[index];
    if (current_chrom.empty())
    {
      current_chrom = ondisc_sptr->getChromatogram(index);
    }

    // fill "dummy" spectrum with chromatogram data
    for (const ChromatogramPeak& cpeak : current_chrom)
    {
      spectrum.emplace_back(cpeak.getRT(), cpeak.getIntensity());
    }

    spectrum.getFloatDataArrays() = current_chrom.getFloatDataArrays();
    spectrum.getIntegerDataArrays() = current_chrom.getIntegerDataArrays();
    spectrum.getStringDataArrays() = current_chrom.getStringDataArrays();

    // Add at least one data point to the chromatogram, otherwise
    // "addPeakLayer" will fail and a segfault occurs later
    if (current_chrom.empty())
    {
      spectrum.emplace_back(-1, 0);
    }
    chrom_exp_sptr->addSpectrum(spectrum);

    // store peptide_sequence if available
    if (current_chrom.getPrecursor().metaValueExists("peptide_sequence"))
    {
      chrom_exp_sptr->setMetaValue("peptide_sequence", current_chrom.getPrecursor().getMetaValue("peptide_sequence"));
    }

    return chrom_exp_sptr;
  }


  Plot1DCanvas::Plot1DCanvas(const Param& preferences, const DIM gravity_axis, QWidget* parent) :
    PlotCanvas(preferences, parent),
    gr_(gravity_axis)
  {
    // for now, default to mz x int
    unit_mapper_ = DimMapper<2>({DIM_UNIT::MZ, DIM_UNIT::INT});

    //Parameter handling
    defaults_.setValue("highlighted_peak_color", "#ff0000", "Highlighted peak color.");
    defaults_.setValue("icon_color", "#000000", "Peak icon color.");
    defaults_.setValue("peak_color", "#0000ff", "Peak color.");
    defaults_.setValue("annotation_color", "#000055", "Annotation color.");
    defaults_.setValue("background_color", "#ffffff", "Background color.");
    defaultsToParam_();
    setName("Plot1DCanvas");
    setParameters(preferences);

    // connect preferences change to the right slot
    connect(this, SIGNAL(preferencesChange()), this, SLOT(currentLayerParamtersChanged_()));
  }

  Plot1DCanvas::~Plot1DCanvas() = default;

  const LayerData1DBase& Plot1DCanvas::getLayer(Size index) const
  {
    return dynamic_cast<const LayerData1DBase&>(layers_.getLayer(index));
  }

  LayerData1DBase& Plot1DCanvas::getLayer(Size index)
  {
    return dynamic_cast<LayerData1DBase&>(layers_.getLayer(index));
  }

  const LayerData1DBase& Plot1DCanvas::getCurrentLayer() const
  {
    return dynamic_cast<const LayerData1DBase&>(layers_.getCurrentLayer());
  }

  LayerData1DBase& Plot1DCanvas::getCurrentLayer()
  {
    return dynamic_cast<LayerData1DBase&>(layers_.getCurrentLayer());
  }

  const DimBase& Plot1DCanvas::getGravityDim() const
  {
    return unit_mapper_.getDim(getGravitator().getGravityAxis());
  }

  const DimBase& Plot1DCanvas::getNonGravityDim() const
  {
    return unit_mapper_.getDim(getGravitator().swap().getGravityAxis());
  }

  bool Plot1DCanvas::addChromLayer(ExperimentSharedPtrType chrom_exp_sptr,
                                   ODExperimentSharedPtrType ondisc_sptr,
                                   OSWDataSharedPtrType chrom_annotation,
                                   const int index,
                                   const String& filename,
                                   const String& basename,
                                   const String& basename_extra)
  {
    // we do not want addChromLayer to trigger repaint, since we have not set the chromatogram data!
    this->blockSignals(true);
    RAIICleanup clean([&]()
    {
      this->blockSignals(false);
    });

    // add chromatogram data as peak spectrum
    if (!PlotCanvas::addChromLayer(std::move(chrom_exp_sptr), std::move(ondisc_sptr), filename))
    {
      return false;
    }
    auto& ld = dynamic_cast<LayerData1DChrom&>(getCurrentLayer());
    ld.setName(basename);
    ld.setNameSuffix(basename_extra);
    ld.getChromatogramAnnotation() = std::move(chrom_annotation); // copy over shared-ptr to OSW-sql data (if available)
    ld.setCurrentIndex(index);                         // use this chrom for visualization
    recalculateRanges_(); // needed here, since 'setCurrentIndex()' changes the current Chromatogram

    setDrawMode(Plot1DCanvas::DM_CONNECTEDLINES);
    //setIntensityMode(Plot1DCanvas::IM_NONE);

    // extend the currently visible area, so the new data is visible
    //auto va = visible_area_.getAreaUnit();

    return true;
  }

  void Plot1DCanvas::activateLayer(Size layer_index)
  {
    layers_.setCurrentLayer(layer_index);

    // no peak is selected
    selected_peak_.clear();

    emit layerActivated(this);
  }

  void Plot1DCanvas::changeVisibleArea1D_(const UnitRange& new_area, bool repaint, bool add_to_stack)
  {
    auto corrected = correctGravityAxisOfVisibleArea_(new_area);
    PlotCanvas::changeVisibleArea_(visible_area_.cloneWith(corrected), repaint, add_to_stack);
  }

  void Plot1DCanvas::changeVisibleArea_(const AreaXYType& new_area, bool repaint, bool add_to_stack)
  {
    changeVisibleArea1D_(visible_area_.cloneWith(new_area).getAreaUnit(), repaint, add_to_stack);
  }
  void Plot1DCanvas::changeVisibleArea_(const UnitRange& new_area, bool repaint, bool add_to_stack)
  {
    changeVisibleArea1D_(new_area, repaint, add_to_stack);
  }
  void Plot1DCanvas::changeVisibleArea_(VisibleArea new_area, bool repaint, bool add_to_stack)
  {
    changeVisibleArea1D_(new_area.getAreaUnit(), repaint, add_to_stack);
  }

  void Plot1DCanvas::dataToWidget(const DPosition<2>& xy_point, QPoint& point, bool flipped)
  {
    Plot1DCanvas::dataToWidget(xy_point.getX(), xy_point.getY(), point, flipped);
  }

  void Plot1DCanvas::dataToWidget(const DPosition<2>& xy_point, DPosition<2>& point, bool flipped)
  {
    QPoint p;
    Plot1DCanvas::dataToWidget(xy_point.getX(), xy_point.getY(), p, flipped);
    point.setX(p.x());
    point.setY(p.y());
  }

  void Plot1DCanvas::dataToWidget(double x, double y, QPoint& point, bool flipped)
  {
    QPoint tmp;
    // adapting gravity dimension is required for percentage mode
    if (gr_.getGravityAxis() == DIM::Y) y *= percentage_factor_;
    else if (gr_.getGravityAxis() == DIM::X) x *= percentage_factor_;
    
    dataToWidget_(x, y, tmp);
    point.setX(tmp.x());
    point.setY(tmp.y());
    
    if (mirror_mode_)
    {
      double alignment_shrink_factor = 1.0;
      if (height() > 10)
      {
        alignment_shrink_factor = (double)(height() - 10) / (double)height();
      }
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
  }

  PointXYType Plot1DCanvas::widgetToData(const QPoint& pos)
  {
    return widgetToData(pos.x(), pos.y());
  }

  PointXYType Plot1DCanvas::widgetToData(double x, double y)
  {
    double actual_y;

    if (mirror_mode_)
    {
      double alignment_shrink_factor = 1.0;
      if (height() > 10)
      {
        alignment_shrink_factor = (double)(height() - 10) / (double)height();
      }
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
    PointXYType p = PlotCanvas::widgetToData_(x, actual_y);
    // adapting gravity dimension is required for percentage mode
    if (gr_.getGravityAxis() == DIM::Y) p.setY(p.getY() / (percentage_factor_));
    else if (gr_.getGravityAxis() == DIM::X) p.setX(p.getX() / (percentage_factor_));

    return p;
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Qt events

  void Plot1DCanvas::mousePressEvent(QMouseEvent* e)
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
          emit sendStatusMessage(String("Measured: d") + getNonGravityDim().getDimNameShort() +  "= " + distance_item->getDistance(), 0);
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
        if (selected_peak_.isValid())
        {
          measurement_start_ = selected_peak_;
          auto peak_xy = getCurrentLayer().peakIndexToXY(measurement_start_, unit_mapper_);
          recalculatePercentageFactor_(getCurrentLayerIndex());
          dataToWidget(peak_xy, measurement_start_point_px_, getCurrentLayer().flipped);
          // use intensity (usually) of mouse, not of the peak
          measurement_start_point_px_ = gr_.gravitateTo(measurement_start_point_px_, last_mouse_pos_);
        }
        else
        {
          measurement_start_.clear();
        }
      }
    }
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Plot1DCanvas::mouseMoveEvent(QMouseEvent* e)
  {
    // mouse position relative to the diagram widget
    QPoint p = e->pos();
    PointXYType data_pos = widgetToData(p);
    emit sendCursorStatus(unit_mapper_.getDim(DIM::X).formattedValue(data_pos[0]),
                          unit_mapper_.getDim(DIM::Y).formattedValue(data_pos[1]));

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
        recalculatePercentageFactor_(getCurrentLayerIndex());
        PointXYType delta = widgetToData(p) - widgetToData(last_mouse_pos_);

        Annotations1DContainer& ann_1d = getCurrentLayer().getCurrentAnnotations();
        for (Annotations1DContainer::Iterator it = ann_1d.begin(); it != ann_1d.end(); ++it)
        {
          if ((*it)->isSelected())
          {
            (*it)->move(delta, gr_, unit_mapper_);         
          }
        }
        update_(OPENMS_PRETTY_FUNCTION);
        last_mouse_pos_ = p;
      }
      else if (action_mode_ == AM_TRANSLATE)
      {
        // translation in data metric
        const double shift = widgetToData(last_mouse_pos_).getX() - widgetToData(p).getX();
        auto new_va = visible_area_.cloneWith(visible_area_.getAreaXY() + PointXYType(shift, 0)).getAreaUnit();

        // change data area
        changeVisibleArea_(new_va);
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
        const auto pixel_area = canvasPixelArea();
        auto r_start = gr_.gravitateMin(last_mouse_pos_, pixel_area);
        auto r_end = gr_.gravitateMax(p, pixel_area);
        rubber_band_.setGeometry(QRect(r_start, r_end).normalized());
        rubber_band_.show(); // if the mouse button is pressed before the zoom key is pressed

        update_(OPENMS_PRETTY_FUNCTION);
      }
    }
    else if (!e->buttons()) // no buttons pressed
    {
      selected_peak_ = near_peak;
      update_(OPENMS_PRETTY_FUNCTION);
    }

    // show coordinates of data arrays
    if (selected_peak_.isValid())
    {
      emit sendStatusMessage(getCurrentLayer().getDataArrayDescription(selected_peak_), 0);
    }
  }

  void Plot1DCanvas::mouseReleaseEvent(QMouseEvent* e)
  {
    if (e->button() == Qt::LeftButton)
    {
      if (action_mode_ == AM_ZOOM)
      {
        rubber_band_.hide();
        QRect rect = rubber_band_.geometry();
        if (rect.width() != 0)
        {
          AreaXYType area(widgetToData(rect.topLeft()), widgetToData(rect.bottomRight()));
          changeVisibleArea_(area, true, true);
        }
      }
      else if (action_mode_ == AM_MEASURE)
      {
        if (selected_peak_.isValid() && measurement_start_.isValid() && selected_peak_.peak != measurement_start_.peak)
        {
          auto start_xy = getCurrentLayer().peakIndexToXY(measurement_start_, unit_mapper_);
          auto end_xy = getCurrentLayer().peakIndexToXY(selected_peak_, unit_mapper_);
          // line should be horizontal at the mouse position --> adapt gravity coordinates
          auto mouse_xy = widgetToData(e->pos());
          start_xy = gr_.gravitateTo(start_xy, mouse_xy);
          end_xy = gr_.gravitateTo(end_xy, mouse_xy);
          recalculatePercentageFactor_(getCurrentLayerIndex());
          // draw line for measured distance between two peaks and annotate with distance in m/z -- use 4 digits to resolve 13C distances between isotopes
          auto* item = new Annotation1DDistanceItem("", start_xy, end_xy);
          item->setText(QString::number(item->getDistance(), 'f', getNonGravityDim().valuePrecision()));
          getCurrentLayer().getCurrentAnnotations().push_front(item);
        }
      }

      moving_annotations_ = false;

      measurement_start_.clear();
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  void Plot1DCanvas::keyPressEvent(QKeyEvent* e)
  {
    // Delete pressed => delete selected annotations from the current layer
    if (e->key() == Qt::Key_Delete)
    {
      e->accept();
      auto peak_layer = dynamic_cast<LayerData1DPeak*>(&getCurrentLayer());
      if (peak_layer)  peak_layer->removePeakAnnotationsFromPeptideHit(getCurrentLayer().getCurrentAnnotations().getSelectedItems());
      getCurrentLayer().getCurrentAnnotations().removeSelectedItems();
      update_(OPENMS_PRETTY_FUNCTION);
    }
    // 'b' pressed && in zoom mode (Ctrl pressed) => select all annotation items
    else if ((e->modifiers() & Qt::ControlModifier) && (e->key() == Qt::Key_B))
    {
      e->accept();
      getCurrentLayer().getCurrentAnnotations().selectAll();
      update_(OPENMS_PRETTY_FUNCTION);
    }
    else
    {
      PlotCanvas::keyPressEvent(e);
    }
  }

  PeakIndex Plot1DCanvas::findPeakAtPosition_(QPoint p)
  {
    //no layers => return invalid peak index
    if (layers_.empty())
    {
      return PeakIndex();
    }
    // mirror mode and p not on same half as active layer => return invalid peak index
    if (mirror_mode_ && (getCurrentLayer().flipped ^ (p.y() > height() / 2)))
    {
      return PeakIndex();
    }
    recalculatePercentageFactor_(getCurrentLayerIndex());

    RangeAllType search_area = unit_mapper_.fromXY(widgetToData(p - QPoint(2, 2)));
    search_area.extend(unit_mapper_.fromXY(widgetToData(p + QPoint(2, 2))));
    return getCurrentLayer().findClosestDataPoint(search_area);
  }

  //////////////////////////////////////////////////////////////////////////////////
  // SLOTS

  void Plot1DCanvas::removeLayer(Size layer_index)
  {
    // remove settings
    layers_.removeLayer(layer_index);
    draw_modes_.erase(draw_modes_.begin() + layer_index);
    peak_penstyle_.erase(peak_penstyle_.begin() + layer_index);

    // update nearest peak
    selected_peak_.clear();

    // abort if there are no layers anymore
    if (layers_.empty())
    {
      overall_data_range_.clearRanges();
      update_(OPENMS_PRETTY_FUNCTION);
      return;
    }

    if (!flippedLayersExist())
    {
      setMirrorModeActive(false);
    }

    // update range area
    recalculateRanges_();
    zoomClear_();
    changeVisibleArea_(overall_data_range_, true, true);
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Plot1DCanvas::setDrawMode(DrawModes mode)
  {
    //no layers
    if (layers_.empty())
    {
      return;
    }
    if (draw_modes_[getCurrentLayerIndex()] != mode)
    {
      draw_modes_[getCurrentLayerIndex()] = mode;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  Plot1DCanvas::DrawModes Plot1DCanvas::getDrawMode() const
  {
    //no layers
    if (layers_.empty())
    {
      return DM_PEAKS;
    }
    return draw_modes_[getCurrentLayerIndex()];
  }

  void Plot1DCanvas::paintEvent(QPaintEvent* e)
  {
    QPainter painter(this);
    paint(&painter, e);
    painter.end();
  }

  void Plot1DCanvas::paint(QPainter* painter, QPaintEvent* e)
  {
    QElapsedTimer timer;
    timer.start();

    // clear
    painter->fillRect(0, 0, this->width(), this->height(),
                      QColor(String(param_.getValue("background_color").toString()).toQString()));

    // we are done if no layer is present
    if (getLayerCount() == 0)
    {
      e->accept();
      return;
    }

    // gridlines
    paintGridLines_(*painter);

    // paint each layer
    for (Size i = 0; i < getLayerCount(); ++i)
    {
      recalculatePercentageFactor_(i);
      auto paint_1d = getLayer(i).getPainter1D();
      paint_1d->paint(painter, this, (int)i);      
    }

    if (show_alignment_)
    {
      drawAlignment_(*painter);
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
        // two x-axes:
        painter->setPen(Qt::black);
        painter->drawLine(0, height() / 2 + 5, width(), height() / 2 + 5);
        painter->drawLine(0, height() / 2 - 5, width(), height() / 2 - 5);
      }
      painter->restore();
    }


    // draw measuring line when in measure mode and valid measurement start peak selected
    if (action_mode_ == AM_MEASURE && measurement_start_.isValid())
    {
      // use start-point + mouse position of non-gravity axis
      QPoint measurement_end_point_px = gr_.swap().gravitateTo(measurement_start_point_px_, last_mouse_pos_);
      auto ps = widgetToData(measurement_start_point_px_);
      auto pe = widgetToData(measurement_end_point_px);
      Annotation1DDistanceItem(QString::number(gr_.swap().gravityDiff(ps, pe), 'f', 4), ps, pe).draw(this, *painter, false);
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

  void Plot1DCanvas::drawHighlightedPeak_(Size layer_index, const PeakIndex& peak, QPainter& painter, bool draw_elongation)
  {
    if (!peak.isValid()) return;
    const auto sel_xy = getLayer(layer_index).peakIndexToXY(peak, unit_mapper_);

    painter.setPen(QPen(QColor(String(param_.getValue("highlighted_peak_color").toString()).toQString()), 2));

    recalculatePercentageFactor_(layer_index);

    QPoint begin;
    dataToWidget(sel_xy, begin, getLayer(layer_index).flipped);

    // paint the cross-hair only for currently selected peaks of the current layer
    if (layer_index == getCurrentLayerIndex() && (peak == measurement_start_ || peak == selected_peak_))
    {
      Painter1DBase::drawCross(begin, &painter, 8);
    }
    // draw elongation as dashed line (while in measure mode and for all existing distance annotations)
    if (draw_elongation)
    {
      QPoint top_end = (getLayer(layer_index).flipped) ? gr_.gravitateMax(begin, canvasPixelArea()) : gr_.gravitateMin(begin, canvasPixelArea());
      Painter1DBase::drawDashedLine(begin, top_end, &painter, String(param_.getValue("highlighted_peak_color").toString()).toQString());
    }
  }

  bool Plot1DCanvas::finishAdding_()
  {
/*    const MSSpectrum& spectrum = getCurrentLayer().getCurrentSpectrum();
    // Abort if no data points are contained (note that all data could be on disk)
    if (spectrum.empty())
    {
      popIncompleteLayer_("Cannot add a dataset that contains no survey scans. Aborting!");
      return false;
    }               
  */

    // add new draw mode and style (default: peaks)
    draw_modes_.push_back(DM_PEAKS);
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
    
    // update nearest peak
    selected_peak_.clear();

    // update ranges
    getCurrentLayer().updateRanges();
    recalculateRanges_();

    resetZoom(false); // no repaint as this is done in setIntensityMode() anyway

    // warn if negative intensities are contained
    if (getCurrentMinIntensity() < 0.0f)
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


  void Plot1DCanvas::drawCoordinates_(QPainter& painter, const PeakIndex& peak)
  {
    if (!peak.isValid())
    {
      return;
    }              
    const auto xy_point = getCurrentLayer().peakIndexToXY(peak, unit_mapper_);
    QStringList lines;
    lines << unit_mapper_.getDim(DIM::X).formattedValue(xy_point.getX()).toQString();
    lines << unit_mapper_.getDim(DIM::Y).formattedValue(xy_point.getY()).toQString();
    drawText_(painter, lines);
  }

  void Plot1DCanvas::drawDeltas_(QPainter& painter, const PeakIndex& start, const PeakIndex& end)
  {
    if (!start.isValid())
    {
      return;
    }

    const auto peak_start = getCurrentLayer().peakIndexToXY(start, unit_mapper_);
    auto peak_end = decltype(peak_start){}; // same as peak_start but without the const
    if (end.isValid())
    {
      peak_end = getCurrentLayer().peakIndexToXY(end, unit_mapper_);
    }
    else
    {
      peak_end = widgetToData_(last_mouse_pos_);
      peak_end = gr_.gravitateNAN(peak_end); // we do not care about the gravity dimension (usually intensity)
    }

    auto dim_text = [](const DimBase& dim, double start_pos, double end_pos, bool ratio /*or difference*/)
    {
      QString result;
      if (ratio)
      {
        result = dim.formattedValue(end_pos / start_pos, " ratio ").toQString();
      }
      else
      {
        result = dim.formattedValue(end_pos - start_pos, " delta ").toQString();
        if (dim.getUnit() == DIM_UNIT::MZ)
        {
          auto ppm = Math::getPPM(end_pos, start_pos);
          result += " (" + QString::number(ppm, 'f', 1) + " ppm)";
        }
      }
      return result;
    };

    QStringList lines;
    lines << dim_text(unit_mapper_.getDim(DIM::X), peak_start.getX(), peak_end.getX(), gr_.getGravityAxis() == DIM::X);
    lines << dim_text(unit_mapper_.getDim(DIM::Y), peak_start.getY(), peak_end.getY(), gr_.getGravityAxis() == DIM::Y);
    drawText_(painter, lines);
  }

  void Plot1DCanvas::recalculatePercentageFactor_(Size layer_index)
  {
    if (intensity_mode_ == IM_PERCENTAGE)
    {
      // maximum value (usually intensity) in whole layer
      const auto max_data_gravity = unit_mapper_.mapRange(getLayer(layer_index).getRange1D()).maxPosition()[(int)gr_.getGravityAxis()];
      percentage_factor_ = 100 / max_data_gravity;
    }
    else
    {
      percentage_factor_ = 1.0;
    }
  }
  void Plot1DCanvas::updateScrollbars_()
  {
    auto xy_overall_area = visible_area_.cloneWith(overall_data_range_1d_).getAreaXY();
    emit updateHScrollbar(xy_overall_area.minPosition()[0], visible_area_.getAreaXY().minPosition()[0], visible_area_.getAreaXY().maxPosition()[0], xy_overall_area.maxPosition()[0]);
    emit updateVScrollbar(1, 1, 1, 1);
  }

  void Plot1DCanvas::horizontalScrollBarChange(int value)
  {
    auto new_area = visible_area_.getAreaXY();
    float shift = value - new_area.center().getX();
    new_area += decltype(new_area)::PositionType(shift, 0);
    changeVisibleArea_(new_area);
  }

  void Plot1DCanvas::showCurrentLayerPreferences()
  {
    Internal::Plot1DPrefDialog dlg(this);
    LayerDataBase& layer = getCurrentLayer();

    ColorSelector* peak_color = dlg.findChild<ColorSelector*>("peak_color");
    ColorSelector* icon_color = dlg.findChild<ColorSelector*>("icon_color");
    ColorSelector* annotation_color = dlg.findChild<ColorSelector*>("annotation_color");
    ColorSelector* bg_color = dlg.findChild<ColorSelector*>("bg_color");
    ColorSelector* selected_color = dlg.findChild<ColorSelector*>("selected_color");

    peak_color->setColor(QColor(String(layer.param.getValue("peak_color").toString()).toQString()));
    icon_color->setColor(QColor(String(layer.param.getValue("icon_color").toString()).toQString()));
    annotation_color->setColor(QColor(String(layer.param.getValue("annotation_color").toString()).toQString()));
    bg_color->setColor(QColor(String(param_.getValue("background_color").toString()).toQString()));
    selected_color->setColor(QColor(String(param_.getValue("highlighted_peak_color").toString()).toQString()));

    if (dlg.exec())
    {
      layer.param.setValue("peak_color", peak_color->getColor().name().toStdString());
      layer.param.setValue("icon_color", icon_color->getColor().name().toStdString());
      layer.param.setValue("annotation_color", annotation_color->getColor().name().toStdString());
      param_.setValue("background_color", bg_color->getColor().name().toStdString());
      param_.setValue("highlighted_peak_color", selected_color->getColor().name().toStdString());

      emit preferencesChange();
    }
  }

  void Plot1DCanvas::currentLayerParamtersChanged_()
  {
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Plot1DCanvas::contextMenuEvent(QContextMenuEvent* e)
  {
    if (layers_.empty()) { return; }

    QMenu* context_menu = new QMenu(this);

    Annotations1DContainer& annots_1d = getCurrentLayer().getCurrentAnnotations();
    Annotation1DItem* annot_item = annots_1d.getItemAt(e->pos());
    bool need_repaint = false; /// will get updated by context menu actions
    if (annot_item)
    {
      annots_1d.deselectAll();
      annots_1d.selectItemAt(e->pos());
      update_(OPENMS_PRETTY_FUNCTION);

      context_menu->addMenu(getCurrentLayer()
        .getContextMenuAnnotation(annot_item, need_repaint));
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
      
      context_menu->addAction("Add label", [&]() {
        addUserLabelAnnotation_(e->pos());
      })->setEnabled(!(mirror_mode_ && (getCurrentLayer().flipped ^ (e->pos().y() > height() / 2))));

      PeakIndex near_peak = findPeakAtPosition_(e->pos());
      context_menu->addAction("Add peak annotation", [&]() {
        addUserPeakAnnotation_(near_peak);
      })->setEnabled(near_peak.isValid());
      
      context_menu->addAction((String("Add peak annotation ") + String(getNonGravityDim().getDimNameShort())).toQString(), [&]() {
        const auto xy_point = getCurrentLayer().peakIndexToXY(near_peak, unit_mapper_);
        QString label = getNonGravityDim().formattedValue(gr_.swap().gravityValue(xy_point)).toQString();
        addPeakAnnotation(near_peak, label, String(getCurrentLayer().param.getValue("peak_color").toString()).toQString());
      })->setEnabled(near_peak.isValid());
      
      context_menu->addSeparator();
      
      context_menu->addAction("Reset alignment", [&]() { 
        resetAlignment();
      })->setEnabled(show_alignment_);

      context_menu->addSeparator();

      context_menu->addAction("Layer meta data", [&]() {
        showMetaData(true);
      });

      QMenu* save_menu = new QMenu("Save");
      
      save_menu->addAction("Layer", [&]() {
        saveCurrentLayer(false);
      });

      save_menu->addAction("Visible layer data", [&]() {
        saveCurrentLayer(true);
      });
      
      save_menu->addAction("As image", [&]() {
        spectrum_widget_->saveAsImage();
      });

      QMenu* settings_menu = new QMenu("Settings");
      
      settings_menu->addAction("Show/hide grid lines", [&]() { 
        showGridLines(!gridLinesShown()); 
      });
      
      settings_menu->addAction("Show/hide axis legends", [&]() {
        emit changeLegendVisibility();
      });
      
      settings_menu->addAction("Style: Stick <--> Area", [&]() {
        if (getDrawMode() != DM_PEAKS)
        {
          setDrawMode(DM_PEAKS);
        }
        else
        {
          setDrawMode(DM_CONNECTEDLINES);
        }
       });

      settings_menu->addAction("Intensity: Absolute <--> Percent", [&]() {
        if (getIntensityMode() != IM_PERCENTAGE)
        {
          setIntensityMode(IM_PERCENTAGE);
        }
        else
        {
          setIntensityMode(IM_SNAP);
        } 
      });

      settings_menu->addAction("Show/hide ion ladder in ID view", [&]() {
        setIonLadderVisible(!isIonLadderVisible());
      });

      settings_menu->addAction("Show/hide automated m/z annotations", [&]() {
        setDrawInterestingMZs(!draw_interesting_MZs_);
      });

      settings_menu->addSeparator();

      settings_menu->addAction("Preferences", [&]() {
        showCurrentLayerPreferences();
      });

      context_menu->addMenu(save_menu);
      context_menu->addMenu(settings_menu);

      // only add to context menu if there is a MS1 map
      auto* peak_layer = dynamic_cast<LayerData1DPeak*>(&getCurrentLayer());
      if (peak_layer)
        {
        if (peak_layer->getPeakData()->containsScanOfLevel(1))
        {
          context_menu->addAction("Switch to 2D view", [&]() {
            emit showCurrentPeaksAs2D();
          });
          context_menu->addAction("Switch to 3D view", [&]() {
            emit showCurrentPeaksAs3D();
          });
        }

        if (peak_layer->getCurrentSpectrum().containsIMData())
        {
          context_menu->addAction("Switch to ion mobility view", [&]() {
            emit showCurrentPeaksAsIonMobility(peak_layer->getCurrentSpectrum());
          });
        }

        if (peak_layer->isDIAData())
        {
          auto l = dynamic_cast<const LayerData1DPeak*>(&getCurrentLayer());
          context_menu->addAction("Switch to DIA-MS view", [&]() {
            emit showCurrentPeaksAsDIA(l->getCurrentSpectrum().getPrecursors()[0], *l->getPeakData().get());
          });
        }
      }

      // add external context menu
      if (context_add_)
      {
        context_menu->addSeparator();
        context_menu->addMenu(context_add_);
      }
    }

    // evaluate menu
    context_menu->exec(mapToGlobal(e->pos()));
    // .. and repaint, depending on action taken
    if (need_repaint)
    {
      update_(OPENMS_PRETTY_FUNCTION);
    }

    e->accept();
  }

  void Plot1DCanvas::setTextBox(const QString& html)
  {
    text_box_content_.setHtml(html);
  }

  void Plot1DCanvas::addUserLabelAnnotation_(const QPoint& screen_position)
  {
    bool ok;
    QString text = QInputDialog::getText(this, "Add label", "Enter text:", QLineEdit::Normal, "", &ok);
    if (ok && !text.isEmpty())
    {
      addLabelAnnotation_(screen_position, text);
    }
  }

  void Plot1DCanvas::addLabelAnnotation_(const QPoint& screen_position, const QString& text)
  {
    recalculatePercentageFactor_(getCurrentLayerIndex());

    PointXYType position = widgetToData(screen_position);
    auto item = new Annotation1DTextItem(position, text);
    getCurrentLayer().getCurrentAnnotations().push_front(item);

    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Plot1DCanvas::addUserPeakAnnotation_(PeakIndex near_peak)
  {
    bool ok;
    QString text = QInputDialog::getText(this, "Add peak annotation", "Enter text:", QLineEdit::Normal, "", &ok);
    if (ok && !text.isEmpty())
    {
      addPeakAnnotation(near_peak, text, QColor(String(getCurrentLayer().param.getValue("peak_color").toString()).toQString()));
    }
  }

  Annotation1DItem* Plot1DCanvas::addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color)
  {
    auto item = getCurrentLayer().addPeakAnnotation(peak_index, text, color);
    update_(OPENMS_PRETTY_FUNCTION);
    return item;
  }

  bool Plot1DCanvas::flippedLayersExist()
  {
    for (Size i = 0; i < getLayerCount(); ++i)
    {
      if (getLayer(i).flipped)
      {
        return true;
      }
    }
    return false;
  }

  void Plot1DCanvas::updateLayer(Size i)
  {
    //update nearest peak
    selected_peak_.clear();

    //update ranges
    recalculateRanges_();

    resetZoom();
    modificationStatus_(i, false);
  }

  void Plot1DCanvas::zoom_(int x, int y, bool zoom_in)
  {
    if (!zoom_in)
    {
      zoomBack_();
    }
    else
    { // only zoom the non-gravity axis
      constexpr PointXYType::CoordinateType zoom_factor = 0.8; // i.e. we crop 20% in total (from left + right, depending on where user clicked)
      // we want to zoom into (x,y), which is in pixel units, hence we need to know the relative position of (x,y) in the widget
      double rel_pos = gr_.getGravityAxis() == DIM::Y ? (PointXYType::CoordinateType)x / width() : (PointXYType::CoordinateType)(height() - y) / height();
      auto new_area = visible_area_.getAreaXY();
      if (gr_.getGravityAxis() == DIM::X) new_area.swapDimensions(); // temporarily swap X<>Y if gravity acts on X
      auto zoomed = Math::zoomIn(new_area.minX(), new_area.maxX(), zoom_factor, rel_pos);
      new_area.setMinX(zoomed.first);
      new_area.setMaxX(zoomed.second);
      if (gr_.getGravityAxis() == DIM::X) new_area.swapDimensions(); // swap back

      if (new_area != visible_area_.getAreaXY())
      {
        zoomAdd_(visible_area_.cloneWith(new_area));
        changeVisibleArea_(*zoom_pos_);  
      }
    }
  }

  /// Go forward in zoom history
  void Plot1DCanvas::zoomForward_()
  {
    // if at end of zoom level then simply add a new zoom
    if (zoom_pos_ == zoom_stack_.end() || (zoom_pos_ + 1) == zoom_stack_.end())
    {
      const auto a = canvasPixelArea();
      zoom_(a.center().getX(), a.center().getY(), true); // calls changeVisibleArea_
      return; 
    }

    // goto next zoom level
    ++zoom_pos_;
    changeVisibleArea_(*zoom_pos_);
  }

  void Plot1DCanvas::translateLeft_(Qt::KeyboardModifiers /*m*/)
  {
    auto xy = visible_area_.getAreaXY();
    // -5% shift in X
    xy -= decltype(xy)::PositionType(0.05 * xy.width(), 0);
    changeVisibleArea_(xy);
  }

  void Plot1DCanvas::translateRight_(Qt::KeyboardModifiers /*m*/)
  {
    auto xy = visible_area_.getAreaXY();
    // +5% shift in X
    xy += decltype(xy)::PositionType(0.05 * xy.width(), 0);
    changeVisibleArea_(xy);
  }

  void Plot1DCanvas::translateForward_()
  {
    auto xy = visible_area_.getAreaXY();
    // +5% shift in Y
    xy += decltype(xy)::PositionType(0, 0.05 * xy.height());
    changeVisibleArea_(xy);
  }

  void Plot1DCanvas::translateBackward_()
  {
    auto xy = visible_area_.getAreaXY();
    // -5% shift in Y
    xy -= decltype(xy)::PositionType(0, 0.05 * xy.height());
    changeVisibleArea_(xy);
  }

  /// Returns whether this widget is currently in mirror mode
  bool Plot1DCanvas::mirrorModeActive() const
  {
    return mirror_mode_;
  }

  /// Sets whether this widget is currently in mirror mode
  void Plot1DCanvas::setMirrorModeActive(bool b)
  {
    mirror_mode_ = b;
    qobject_cast<Plot1DWidget*>(spectrum_widget_)->toggleMirrorView(b);
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Plot1DCanvas::paintGridLines_(QPainter& painter)
  {
    if (!show_grid_ || !spectrum_widget_)
    {
      return;
    }
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

      for (const auto& line : spectrum_widget_->xAxis()->gridLines()[j])
      {
        int x = static_cast<int>(Math::intervalTransformation(line, spectrum_widget_->xAxis()->getAxisMinimum(), spectrum_widget_->xAxis()->getAxisMaximum(), xl, xh));
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

      for (const auto& line : spectrum_widget_->yAxis()->gridLines()[j])
      {
        int y = static_cast<int>(Math::intervalTransformation(line, spectrum_widget_->yAxis()->getAxisMinimum(), spectrum_widget_->yAxis()->getAxisMaximum(), yl, yh));
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

  void Plot1DCanvas::performAlignment(Size layer_index_1, Size layer_index_2, const Param& param)
  {
    alignment_layer_1_ = layer_index_1;
    alignment_layer_2_ = layer_index_2;
    aligned_peaks_mz_delta_.clear();
    aligned_peaks_indices_.clear();

    if (layer_index_1 >= getLayerCount() || layer_index_2 >= getLayerCount())
    {
      return;
    }
    auto ptr_layer_1 = dynamic_cast<const LayerData1DPeak*>(&getLayer(layer_index_1));
    auto ptr_layer_2 = dynamic_cast<const LayerData1DPeak*>(&getLayer(layer_index_2));
    if (ptr_layer_1 == nullptr || ptr_layer_2 == nullptr)
    {
      return;
    }
    const ExperimentType::SpectrumType& spectrum_1 = ptr_layer_1->getCurrentSpectrum();
    const ExperimentType::SpectrumType& spectrum_2 = ptr_layer_2->getCurrentSpectrum();

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

  void Plot1DCanvas::resetAlignment()
  {
    aligned_peaks_indices_.clear();
    aligned_peaks_mz_delta_.clear();
    qobject_cast<Plot1DWidget*>(spectrum_widget_)->resetAlignment();
    show_alignment_ = false;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Plot1DCanvas::drawAlignment_(QPainter& painter)
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
      auto ptr_layer_1 = dynamic_cast<const LayerData1DPeak*>(&getLayer(alignment_layer_1_));
      if (ptr_layer_1 == nullptr)
      {
        return;
      }
      const ExperimentType::SpectrumType& spectrum_1 = ptr_layer_1->getCurrentSpectrum();
      recalculatePercentageFactor_(alignment_layer_1_);
      for (Size i = 0; i < getAlignmentSize(); ++i)
      {
        dataToWidget(spectrum_1[aligned_peaks_indices_[i].first].getMZ(), 0, begin_p, false);
        dataToWidget(spectrum_1[aligned_peaks_indices_[i].first].getMZ(), spectrum_1[aligned_peaks_indices_[i].first].getIntensity(), end_p, false);
        painter.drawLine(begin_p.x(), begin_p.y(), end_p.x(), end_p.y());
      }
    }
    painter.restore();
  }

  Size Plot1DCanvas::getAlignmentSize()
  {
    return aligned_peaks_mz_delta_.size();
  }

  double Plot1DCanvas::getAlignmentScore() const
  {
    return alignment_score_;
  }

  void Plot1DCanvas::intensityModeChange_()
  {
    changeVisibleArea_(visible_area_, false, false); // updates y-axis
    ensureAnnotationsWithinDataRange_();
    // update axes (e.g. make it Log-scale)
    if (spectrum_widget_)
    {
      spectrum_widget_->updateAxes();
    }
    update_(OPENMS_PRETTY_FUNCTION);
  }

  RangeAllType Plot1DCanvas::correctGravityAxisOfVisibleArea_(UnitRange area)
  {
    // depending on intensity mode, the y-axis either shows
    //  The maximum range (normal & log) 
    //  or the local maximum (snap mode)
    //  or [0, 100] (percentage mode)
    if (intensity_mode_ == PlotCanvas::IntensityModes::IM_SNAP)
    {                                        // find the range of the current data (as determined by x-axis)
      area.clear(getGravityDim().getUnit()); // delete gravity (e.g. intensity), only keep the non-gravity range (e.g. m/z)
      for (Size i = 0; i < getLayerCount(); ++i)
      {
        area.extend(getLayer(i).getRangeForArea(area));
      }
      auto& intensity = getGravityDim().map(area); // make sure y-axis spans [0, max * TOP_MARGIN]
      intensity.setMin(0); // make sure we start at 0
      intensity.extend(intensity.getMax() * TOP_MARGIN);
    }
    else if (intensity_mode_ == PlotCanvas::IntensityModes::IM_PERCENTAGE)
    {
      auto& intensity = getGravityDim().map(area);
      intensity = RangeBase(0, 100 * TOP_MARGIN);
    }
    else
    { // use y-range of all layers
      auto& intensity = getGravityDim().map(area);
      intensity = getGravityDim().map(overall_data_range_1d_);
      intensity.setMin(0); // make sure we start at 0
    }
    return area;
  }

  void Plot1DCanvas::ensureAnnotationsWithinDataRange_()
  {
    for (Size i = 0; i < getLayerCount(); ++i)
    {
      recalculatePercentageFactor_(i);
      Annotations1DContainer& ann_1d = getLayer(i).getCurrentAnnotations();
      for (Annotations1DContainer::Iterator it = ann_1d.begin(); it != ann_1d.end(); ++it)
      {
        (*it)->ensureWithinDataRange(this, i);
      }
    }
  }

  void Plot1DCanvas::flipLayer(Size index)
  {
    if (index < getLayerCount())
    {
      getLayer(index).flipped = !getLayer(index).flipped;
    }
  }

  void Plot1DCanvas::activateSpectrum(Size index, bool repaint)
  {
    // clear selected peak, so we do not accidentally access an invalid index next time when moving the mouse
    selected_peak_.clear();
    
    if (getCurrentLayer().hasIndex(index))
    {
      getCurrentLayer().setCurrentIndex(index);
      recalculateRanges_(); // adapt overall_range_(1d)_
      changeVisibleArea_(visible_area_, repaint, false); // updates y-axis based on new spectrum
    }
  }

  void Plot1DCanvas::setCurrentLayerPeakPenStyle(Qt::PenStyle ps)
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

  std::vector<std::pair<Size, Size> > Plot1DCanvas::getAlignedPeaksIndices()
  {
    return aligned_peaks_indices_;
  }

  void Plot1DCanvas::setIonLadderVisible(bool show)
  {
    if (ion_ladder_visible_ != show)
    {
      ion_ladder_visible_ = show;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  void Plot1DCanvas::setDrawInterestingMZs(bool enable)
  {
    if (draw_interesting_MZs_ != enable)
    {
      draw_interesting_MZs_ = enable;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  bool Plot1DCanvas::isIonLadderVisible() const
  {
    return ion_ladder_visible_;
  }

  bool Plot1DCanvas::isDrawInterestingMZs() const
  {
    return draw_interesting_MZs_;
  }

} //Namespace
