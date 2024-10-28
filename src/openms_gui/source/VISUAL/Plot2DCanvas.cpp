// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/DIALOGS/FeatureEditDialog.h>
#include <OpenMS/VISUAL/DIALOGS/Plot2DPrefDialog.h>
#include <OpenMS/VISUAL/INTERFACES/IPeptideIds.h>
#include <OpenMS/VISUAL/LayerDataConsensus.h>
#include <OpenMS/VISUAL/LayerDataFeature.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/Painter2DBase.h>
#include <OpenMS/VISUAL/Plot2DCanvas.h>
#include <OpenMS/VISUAL/PlotWidget.h>
//STL
#include <algorithm>

//QT
#include <QMouseEvent>
#include <QPainter>
#include <QPolygon>
#include <QElapsedTimer>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMessageBox>

#define PEN_SIZE_MAX_LIMIT 100    // maximum size of a rectangle representing a point for raw peak data
#define PEN_SIZE_MIN_LIMIT 1      // minimum. This should not be changed without adapting the way dots are plotted
                                  // (might lead to inconsistencies when switching between drawing modes of
                                  //  paintMaximumIntensities() vs. paintAllIntensities() )

// the two following constants describe the valid range of the 'canvas_coverage_min_' member (adaptable by the user)
#define CANVAS_COVERAGE_MIN_LIMITHIGH 0.5
#define CANVAS_COVERAGE_MIN_LIMITLOW 0.1


using namespace std;

namespace OpenMS
{
  using namespace Internal;

  Plot2DCanvas::Plot2DCanvas(const Param & preferences, QWidget * parent) :
    PlotCanvas(preferences, parent),
    pen_size_min_(1),
    pen_size_max_(20),
    canvas_coverage_min_(0.2)
  {
    //Parameter handling
    defaults_.setValue("background_color", "#ffffff", "Background color.");
    defaults_.setValue("interpolation_steps", 1000, "Number of interpolation steps for peak gradient pre-calculation.");
    defaults_.setMinInt("interpolation_steps", 1);
    defaults_.setMaxInt("interpolation_steps", 1000);
    defaults_.setValue("dot:gradient", "Linear|0,#eeeeee;1,#ffea00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000", "Multi-color gradient for peaks.");
    defaults_.setValue("dot:feature_icon", "circle", "Icon used for features and consensus features.");
    defaults_.setValidStrings("dot:feature_icon", {"diamond","square","circle","triangle"});
    defaults_.setValue("dot:feature_icon_size", 4, "Icon size used for features and consensus features.");
    defaults_.setMinInt("dot:feature_icon_size", 1);
    defaults_.setMaxInt("dot:feature_icon_size", 999);
    defaultsToParam_();
    setName("Plot2DCanvas");
    setParameters(preferences);

    linear_gradient_.fromString(param_.getValue("dot:gradient"));

    // connect preferences change to the right slot
    connect(this, SIGNAL(preferencesChange()), this, SLOT(currentLayerParametersChanged_()));
  }

  Plot2DCanvas::~Plot2DCanvas() = default;

  void Plot2DCanvas::highlightPeak_(QPainter & painter, const PeakIndex & peak)
  {
    if (!peak.isValid())
    {
      return;
    }
    //determine coordinates;
    auto pos_xy = getCurrentLayer().peakIndexToXY(peak, unit_mapper_);

    // paint highlighted peak
    painter.save();
    painter.setPen(QPen(Qt::red, 2));

    auto p_px = dataToWidget_(pos_xy);
    painter.drawEllipse(p_px.x() - 5, p_px.y() - 5, 10, 10);

    //restore painter
    painter.restore();
  }

  PeakIndex Plot2DCanvas::findNearestPeak_(const QPoint& pos)
  {
    // no layers => return invalid peak index
    if (layers_.empty())
    {
      return PeakIndex();
    }
    // constructing area around the current mouse position
    auto a = visible_area_;
    a.setArea(VisibleArea::AreaXYType(widgetToData_(pos - QPoint(5, 5)), widgetToData_(pos + QPoint(5, 5))));
    return getCurrentLayer().findHighestDataPoint(a.getAreaUnit());
  }
  
  double Plot2DCanvas::adaptPenScaling_(double ratio_data2pixel, double& pen_width) const
  {
    // is the coverage OK using current pen width?
    bool has_low_pixel_coverage_withpen = ratio_data2pixel*pen_width < canvas_coverage_min_;
    int merge_factor(1);
    if (has_low_pixel_coverage_withpen)
    { // scale up the sparse dimension until we reach the desired coverage (this will lead to overlap in the crowded dimension)
      double scale_factor = canvas_coverage_min_ / ratio_data2pixel;
      // however, within bounds (no need to check the pen_size_min_, because we can only exceed here, not underestimate)
      scale_factor = std::min(pen_size_max_, scale_factor);
      // The difference between the original pen_width vs. this scale
      // gives the number of peaks to merge in the crowded dimension
      merge_factor = scale_factor / pen_width;
      // set pen width to the new scale
      pen_width = scale_factor;
    }
    return merge_factor;
  }


  void Plot2DCanvas::intensityModeChange_()
  {
    String gradient_str;
    if (intensity_mode_ == IM_LOG)
    {
      gradient_str = MultiGradient::getDefaultGradientLogarithmicIntensityMode().toString();
    }
    else // linear
    {
      gradient_str = linear_gradient_.toString();
    }
    if (layers_.empty())
    {
      return;
    }
    layers_.getCurrentLayer().param.setValue("dot:gradient", gradient_str);
    for (Size i = 0; i < layers_.getLayerCount(); ++i)
    {
      recalculateDotGradient_(i);
    }

    PlotCanvas::intensityModeChange_();
  }

  void Plot2DCanvas::recalculateDotGradient_(Size layer)
  {
    getLayer(layer).gradient.fromString(getLayer(layer).param.getValue("dot:gradient"));
    if (intensity_mode_ == IM_LOG)
    {
      getLayer(layer).gradient.activatePrecalculationMode(0.0, std::log1p(overall_data_range_.getMaxIntensity()), param_.getValue("interpolation_steps"));
    }
    else
    {
      getLayer(layer).gradient.activatePrecalculationMode(0.0, overall_data_range_.getMaxIntensity(), param_.getValue("interpolation_steps"));
    }
  }

  void Plot2DCanvas::recalculateCurrentLayerDotGradient()
  {
    recalculateDotGradient_(layers_.getCurrentLayerIndex());
  }

  void Plot2DCanvas::pickProjectionLayer()
  {
    // find the last (visible) peak layers
    Size layer_count = 0;
    Size last_layer = 0;
    Size visible_layer_count = 0;
    Size visible_last_layer = 0;
    for (Size i = 0; i < getLayerCount(); ++i)
    {
      if (getLayer(i).type == LayerDataBase::DT_PEAK)
      {
        layer_count++;
        last_layer = i;

        if (getLayer(i).visible)
        {
          visible_layer_count++;
          visible_last_layer = i;
        }
      }
      if (getLayer(i).type == LayerDataBase::DT_CHROMATOGRAM)
      {
        //TODO CHROM
      }
    }

    // try to find the right layer to project
    const LayerDataBase* layer = nullptr;
    //first choice: current layer
    if (layer_count != 0 && getCurrentLayer().type == LayerDataBase::DT_PEAK)
    {
      layer = &(getCurrentLayer());
    }
    //second choice: the only peak layer
    else if (layer_count == 1)
    {
      layer = &(getLayer(last_layer));
    }
    //third choice: the only visible peak layer
    else if (visible_layer_count == 1)
    {
      layer = &(getLayer(visible_last_layer));
    }
    //no layer with peaks: disable projections
    else
    {
      emit toggleProjections();
      return;
    }

    emit showProjections(layer);    
  }

  bool Plot2DCanvas::finishAdding_()
  {
    // deselect all peaks
    selected_peak_.clear();
    measurement_start_.clear();

    auto& layer = getCurrentLayer();
    layer.updateRanges(); // required for minIntensity() below and hasRange()
    if (layer.getRange().hasRange() == HasRangeType::NONE)
    {
      popIncompleteLayer_("Cannot add a dataset that contains no data. Aborting!");
      return false;
    }
    // overall values update
    recalculateRanges_();

    // pick dimensions to show (based on data)


    update_buffer_ = true;

    if (getLayerCount() == 1)
    {
      resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway
    }
    else if (getLayerCount() == 2)
    {
      setIntensityMode(IM_PERCENTAGE);
    }
    intensityModeChange_();

    emit layerActivated(this);

    // warn if negative intensities are contained
    if (getCurrentMinIntensity() < 0)
    {
      QMessageBox::warning(this, "Warning", "This dataset contains negative intensities. Use it at your own risk!");
    }
    
    return true;
  }

  void Plot2DCanvas::removeLayer(Size layer_index)
  {
    if (layer_index >= getLayerCount())
    {
      return;
    }

    // remove the data
    layers_.removeLayer(layer_index);

    // update visible area and boundaries
    auto old_data_range = overall_data_range_;
    recalculateRanges_();

    // only reset zoom if data range has been changed
    if (old_data_range != overall_data_range_)
    {
      resetZoom(false); // no repaint as this is done in intensityModeChange_() anyway
    }

    if (layers_.empty())
    {
      overall_data_range_.clearRanges();
      update_buffer_ = true;
      update_(OPENMS_PRETTY_FUNCTION);
      return;
    }

    // unselect all peaks
    selected_peak_.clear();
    measurement_start_.clear();

    intensityModeChange_();

    emit layerActivated(this);
  }

  // change the current layer
  void Plot2DCanvas::activateLayer(Size layer_index)
  {
    // unselect all peaks
    selected_peak_.clear();
    measurement_start_.clear();

    layers_.setCurrentLayer(layer_index);
    emit layerActivated(this);

    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Plot2DCanvas::recalculateSnapFactor_()
  {
    snap_factors_ = vector<double>(getLayerCount(), 1.0);

    if (intensity_mode_ == IM_SNAP)
    {
      for (Size i = 0; i < getLayerCount(); i++)
      {
        if (getLayer(i).visible)
        {
          auto local_max  = -numeric_limits<Peak1D::IntensityType>::max();
          if (auto* lp = dynamic_cast<LayerDataPeak*>(&getLayer(i)))
          {
            for (ExperimentType::ConstAreaIterator it = lp->getPeakData()->areaBeginConst(visible_area_.getAreaUnit().getMinRT(), visible_area_.getAreaUnit().getMaxRT(),
                                                                                          visible_area_.getAreaUnit().getMinMZ(), visible_area_.getAreaUnit().getMaxMZ());
                 it != lp->getPeakData()->areaEndConst();
                 ++it)
            {
              PeakIndex pi = it.getPeakIndex();
              if (it->getIntensity() > local_max && getLayer(i).filters.passes((*lp->getPeakData())[pi.spectrum], pi.peak))
              {
                local_max = it->getIntensity();
              }
            }
          }
          else if (auto* lp = dynamic_cast<LayerDataFeature*>(&getLayer(i))) // features
          {
            for (FeatureMapType::ConstIterator it = lp->getFeatureMap()->begin(); it != lp->getFeatureMap()->end();
                 ++it)
            {
              if (visible_area_.getAreaUnit().containsRT(it->getRT()) &&
                  visible_area_.getAreaUnit().containsMZ(it->getMZ()) &&
                  getLayer(i).filters.passes(*it))
              {
                local_max = std::max(local_max, it->getIntensity());
              }
            }
          }
          else if (auto* lp = dynamic_cast<LayerDataConsensus*>(&getLayer(i))) // consensus
          {
            for (ConsensusMapType::ConstIterator it = lp->getConsensusMap()->begin(); it != lp->getConsensusMap()->end();
                 ++it)
            {
              if (visible_area_.getAreaUnit().containsRT(it->getRT()) &&
                  visible_area_.getAreaUnit().containsMZ(it->getMZ()) &&
                  getLayer(i).filters.passes(*it) &&
                  it->getIntensity() > local_max)
              {
                local_max = it->getIntensity();
              }
            }
          }
          else if (getLayer(i).type == LayerDataBase::DT_CHROMATOGRAM)  // chromatogram
          {
            //TODO CHROM
          }
          else if (getLayer(i).type == LayerDataBase::DT_IDENT)         // identifications
          {
            //TODO IDENT
          }

          if (local_max > 0.0)
          {
            snap_factors_[i] = overall_data_range_.getMaxIntensity() / local_max;
          }
        }
      }
    }
  }

  void Plot2DCanvas::updateScrollbars_()
  {
    GenericArea ga(&unit_mapper_);
    auto all = ga.setArea(overall_data_range_).getAreaXY();
    const auto& vis = visible_area_.getAreaXY();
    emit updateHScrollbar(all.minX(), vis.minX(), vis.maxX(), all.maxX());
    emit updateVScrollbar(all.minY(), vis.minY(), vis.maxY(), all.maxY());
  }

  void Plot2DCanvas::horizontalScrollBarChange(int value)
  {
    auto new_area = visible_area_;
    auto new_XY = new_area.getAreaXY();
    auto X_width = new_XY.width();
    new_XY.setMinX(value);
    new_XY.setMaxX(value + X_width);
    new_area.setArea(new_XY);
    changeVisibleArea_(new_area);
  }

  void Plot2DCanvas::verticalScrollBarChange(int value)
  {
    // invert 'value' (since the VERTICAL(!) scrollbar's range is negative -- see PlotWidget::updateVScrollbar())
    value *= -1;
    auto new_area = visible_area_;
    auto new_XY = new_area.getAreaXY();
    auto Y_height = new_XY.height();
    new_XY.setMinY(value);
    new_XY.setMaxY(value + Y_height);
    new_area.setArea(new_XY);
    changeVisibleArea_(new_area);
  }

  void Plot2DCanvas::paintEvent(QPaintEvent * e)
  {
    //Only fill background if no layer is present
    if (getLayerCount() == 0)
    {
      QPainter painter;
      painter.begin(this);
      painter.fillRect(0, 0, this->width(), this->height(), QColor(String(param_.getValue("background_color").toString()).toQString()));
      painter.end();
      e->accept();
      return;
    }

#ifdef DEBUG_TOPPVIEW
    cout << "BEGIN " << OPENMS_PRETTY_FUNCTION << endl;
    cout << "  Visible area -- m/z: " << visible_area_.minX() << " - " << visible_area_.maxX() << " rt: " << visible_area_.minY() << " - " << visible_area_.maxY() << endl;
    cout << "  Overall area -- m/z: " << overall_data_range_.minPosition()[0] << " - " << overall_data_range_.maxPosition()[0] << " rt: " << overall_data_range_.minPosition()[1] << " - " << overall_data_range_.maxPosition()[1] << endl;
#endif

    //timing
    QElapsedTimer overall_timer;
    if (show_timing_)
    {
      overall_timer.start();
      if (update_buffer_)
      {
        cout << "Updating buffer:" << endl;
      }
      else
      {
        cout << "Copying buffer:" << endl;
      }
    }

    QPainter painter;
    if (update_buffer_)
    {
      update_buffer_ = false;

      // recalculate snap factor
      recalculateSnapFactor_();

      buffer_.fill(QColor(String(param_.getValue("background_color").toString()).toQString()).rgb());
      painter.begin(&buffer_);
      QElapsedTimer layer_timer;

      for (Size i = 0; i < getLayerCount(); i++)
      {
        // timing
        if (show_timing_)
        {
          layer_timer.start();
        }
        if (getLayer(i).visible)
        {
          // update factors (snap and percentage)
          percentage_factor_ = 1.0;
          if (intensity_mode_ == IM_PERCENTAGE)
          {
            if (getLayer(i).getMaxIntensity() > 0)
            {
              percentage_factor_ = overall_data_range_.getMaxIntensity() / getLayer(i).getMaxIntensity();
            }
          }
          getLayer(i).getPainter2D()->paint(&painter, this, i);
        }
        // timing
        if (show_timing_)
        {
          cout << "  -layer " << i << " time: " << layer_timer.elapsed() << " ms" << endl;
        }
      }
      paintGridLines_(painter);
      painter.end();
    }

    painter.begin(this);

    // copy peak data from buffer
     /*
         * Suppressed warning QVector<QRect> QRegion::rects() const is deprecated
         * Use begin()/end() instead, from Qt 5.8
         */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    QVector<QRect> rects = e->region().rects();
#pragma GCC diagnostic pop
    for (int i = 0; i < (int)rects.size(); ++i)
    {
      painter.drawImage(rects[i].topLeft(), buffer_, rects[i]);
    }

    // draw measurement peak
    if (action_mode_ == AM_MEASURE && measurement_start_.isValid())
    {
      painter.setPen(Qt::black);

      QPoint line_begin;
      // start of line
      if (selected_peak_.isValid())
      {
        auto data_xy = getCurrentLayer().peakIndexToXY(selected_peak_, unit_mapper_);
        line_begin = dataToWidget_(data_xy);
      }
      else
      {
        line_begin = last_mouse_pos_;
      }

      // end of line
      auto data_xy = getCurrentLayer().peakIndexToXY(measurement_start_, unit_mapper_);
      auto line_end = dataToWidget_(data_xy);
      painter.drawLine(line_begin, line_end);

      highlightPeak_(painter, measurement_start_);
    }

    // draw convex hulls or consensus feature elements
    if (selected_peak_.isValid())
    {
      getCurrentLayer().getPainter2D()->highlightElement(&painter, this, selected_peak_);
    }

    if (action_mode_ == AM_MEASURE || action_mode_ == AM_TRANSLATE)
    {
      highlightPeak_(painter, selected_peak_);
    }

    // draw delta for measuring
    if (action_mode_ == AM_MEASURE && measurement_start_.isValid())
    {
      drawDeltas_(painter, measurement_start_, selected_peak_);
    }
    else
    {
      drawCoordinates_(painter, selected_peak_);
    }

    painter.end();
    if (show_timing_)
    {
      cout << "  -overall time: " << overall_timer.elapsed() << " ms" << endl << endl;
    }
  }

  void Plot2DCanvas::drawCoordinates_(QPainter & painter, const PeakIndex & peak)
  {
    if (!peak.isValid()) return;

    const auto xy_point = getCurrentLayer().peakIndexToXY(peak, unit_mapper_);
    QStringList lines;
    lines << unit_mapper_.getDim(DIM::X).formattedValue(xy_point.getX()).toQString();
    lines << unit_mapper_.getDim(DIM::Y).formattedValue(xy_point.getY()).toQString();
    if (unit_mapper_.getDim(DIM::X).getUnit() != DIM_UNIT::INT && unit_mapper_.getDim(DIM::Y).getUnit() != DIM_UNIT::INT)
    { // if intensity is not mapped to X or Y, add it
      // Note: it may be cleaner to hoist this function into the derived classes of Painter2D, 
      //       if the logic here depends on the actual Layer type (currently, 'INT' should work fine for all).
      DimMapper<2> int_mapper({DIM_UNIT::INT, DIM_UNIT::INT});
      const auto int_point = getCurrentLayer().peakIndexToXY(peak, int_mapper);
      lines << int_mapper.getDim(DIM::X).formattedValue(int_point.getX()).toQString();
    }
    drawText_(painter, lines);
  }

  void Plot2DCanvas::drawDeltas_(QPainter & painter, const PeakIndex & start, const PeakIndex & end)
  {
    if (!start.isValid())
    {
      return;
    }

    // mapper obtain intensity from a PeakIndex
    DimMapper<2> intensity_mapper({DIM_UNIT::INT, DIM_UNIT::INT});


    const auto peak_start = getCurrentLayer().peakIndexToXY(start, unit_mapper_);
    const auto peak_start_int = getCurrentLayer().peakIndexToXY(start, intensity_mapper);
    auto peak_end = decltype(peak_start) {}; // same as peak_start but without the const
    auto peak_end_int = decltype(peak_start) {}; // same as peak_start but without the const
    if (end.isValid())
    {
      peak_end = getCurrentLayer().peakIndexToXY(end, unit_mapper_);
      peak_end_int = getCurrentLayer().peakIndexToXY(end, intensity_mapper);
    }
    else
    {
      peak_end = widgetToData_(last_mouse_pos_);
      peak_end_int = decltype(peak_end_int) {};
    }

    auto dim_text = [](const DimBase& dim, double start_pos, double end_pos, bool ratio /*or difference*/) {
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
    lines << dim_text(unit_mapper_.getDim(DIM::X), peak_start.getX(), peak_end.getX(), false);
    lines << dim_text(unit_mapper_.getDim(DIM::Y), peak_start.getY(), peak_end.getY(), false);
    lines << dim_text(DimINT(), peak_start_int.getY(), peak_end_int.getY(), true); // ratio
    drawText_(painter, lines);
  }

  void Plot2DCanvas::mousePressEvent(QMouseEvent * e)
  {
    last_mouse_pos_ = e->pos();

    if (e->button() == Qt::LeftButton)
    {
      if (action_mode_ == AM_MEASURE)
      {
        if (selected_peak_.isValid())
        {
          measurement_start_ = selected_peak_;
        }
        else
        {
          measurement_start_.clear();
        }
      }
      else if (action_mode_ == AM_ZOOM)
      {
        //translate (if not moving features)
        if ( !(getCurrentLayer().type == LayerDataBase::DT_FEATURE) || !selected_peak_.isValid())
        {
          rubber_band_.setGeometry(QRect(e->pos(), QSize()));
          rubber_band_.show();
        }
      }
    }
  }

  void Plot2DCanvas::mouseMoveEvent(QMouseEvent* e)
  {
    grabKeyboard();     // (re-)grab keyboard after it has been released by unhandled key
    QPoint pos = e->pos();
    PointXYType data_pos = widgetToData_(pos);
    emit sendCursorStatus(unit_mapper_.getDim(DIM::X).formattedValue(data_pos[0]),
                          unit_mapper_.getDim(DIM::Y).formattedValue(data_pos[1]));

    PeakIndex near_peak = findNearestPeak_(pos);

    //highlight current peak and display peak coordinates
    if (action_mode_ == AM_MEASURE || (action_mode_ == AM_TRANSLATE && !(e->buttons() & Qt::LeftButton)))
    {
      //highlight peak
      selected_peak_ = near_peak;
      update_(OPENMS_PRETTY_FUNCTION);

      //show meta data in status bar (if available)
      if (selected_peak_.isValid())
      {
        String status;
        auto* lf = dynamic_cast<LayerDataFeature*>(&getCurrentLayer());
        auto* lc = dynamic_cast<LayerDataConsensus*>(&getCurrentLayer());
        if (lf || lc)
        {
          //add meta info
          const BaseFeature* f;
          if (lf)
          {
            f = &selected_peak_.getFeature(*lf->getFeatureMap());
          }
          else
          {
            f = &selected_peak_.getFeature(*lc->getConsensusMap());
          }
          std::vector<String> keys;
          f->getKeys(keys);
          for (Size m = 0; m < keys.size(); ++m)
          {
            status += " " + keys[m] + ": ";
            const DataValue& dv = f->getMetaValue(keys[m]);
            if (dv.valueType() == DataValue::DOUBLE_VALUE)
            { // use less precision for large numbers, for better readability
              int precision(2);
              if ((double)dv < 10) precision = 5;
              // ... and add 1k separators, e.g. '540,321.99'
              status += QLocale::c().toString((double)dv, 'f', precision);
            }
            else
            {
              status += (String)dv;
            }
          }
        }
        else if (auto* lp = dynamic_cast<LayerDataPeak*>(&getCurrentLayer()))
        {
          //meta info
          const ExperimentType::SpectrumType & s = selected_peak_.getSpectrum(*lp->getPeakData());
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
        }
        else if (getCurrentLayer().type == LayerDataBase::DT_CHROMATOGRAM) // chromatogram
        {
          //TODO CHROM
        }
        if (status != "")
        {
          emit sendStatusMessage(status, 0);
        }
      }
    }
    else if (action_mode_ == AM_ZOOM)
    {
      //Zoom mode => no peak should be selected
      selected_peak_.clear();
      update_(OPENMS_PRETTY_FUNCTION);
    }

    if (action_mode_ == AM_MEASURE)
    {
      last_mouse_pos_ = pos;
    }
    else if (action_mode_ == AM_ZOOM)
    {
      //if mouse button is held down, enlarge the selection
      if (e->buttons() & Qt::LeftButton)
      {
        rubber_band_.setGeometry(QRect(last_mouse_pos_, pos).normalized());
        rubber_band_.show();         //if the mouse button is pressed before the zoom key is pressed

        update_(OPENMS_PRETTY_FUNCTION);
      }
    }
    else if (action_mode_ == AM_TRANSLATE)
    {
      if (e->buttons() & Qt::LeftButton)
      {
        // move feature
        auto* lf = dynamic_cast<LayerDataFeature*>(&getCurrentLayer());
        if (getCurrentLayer().modifiable && lf && selected_peak_.isValid())
        {
          RangeType dr;
          unit_mapper_.fromXY(widgetToData_(pos), dr);
          // restrict the movement to the data range
          dr.pushInto(overall_data_range_);
          (*lf->getFeatureMap())[selected_peak_.peak].setRT(dr.getMinRT());
          (*lf->getFeatureMap())[selected_peak_.peak].setMZ(dr.getMinMZ());

          update_buffer_ = true;
          update_(OPENMS_PRETTY_FUNCTION);
          modificationStatus_(layers_.getCurrentLayerIndex(), true);
        }
        else // translate
        { // calculate data coordinates of shift
          PointXYType old_data = widgetToData_(last_mouse_pos_);
          PointXYType new_data = widgetToData_(pos);
          // compute new area
          auto new_visible_area = visible_area_;
          new_visible_area.setArea(new_visible_area.getAreaXY() + (old_data - new_data));
          // publish (bounds checking is done inside changeVisibleArea_)
          changeVisibleArea_(new_visible_area);
          last_mouse_pos_ = pos;
        }
      }
    }
  }

  void Plot2DCanvas::mouseReleaseEvent(QMouseEvent * e)
  {
    if (e->button() == Qt::LeftButton)
    {
      if (action_mode_ == AM_MEASURE)
      {
        if (!selected_peak_.isValid())
        {
          measurement_start_.clear();
        }
        measurement_start_.clear();
        update_(OPENMS_PRETTY_FUNCTION);
      }
      else if (action_mode_ == AM_ZOOM)
      {
        rubber_band_.hide();
        QRect rect = rubber_band_.geometry();
        if (rect.width() != 0 && rect.height() != 0)
        {
          auto new_area = visible_area_.cloneWith(AreaXYType(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight())));
          changeVisibleArea_(new_area, true, true);
        }
      }
    }
  }

  void Plot2DCanvas::contextMenuEvent(QContextMenuEvent * e)
  {
    // abort if there are no layers
    if (layers_.empty())
    {
      return;
    }
    const LayerDataBase& layer = getCurrentLayer();

    QMenu* context_menu = new QMenu(this);

    QAction* a = nullptr;
    QAction* result = nullptr;

    //Display name and warn if current layer invisible
    String layer_name = String("Layer: ") + layer.getName();
    if (!layer.visible)
    {
      layer_name += " (invisible)";
    }
    context_menu->addAction(layer_name.toQString())->setEnabled(false);
    context_menu->addSeparator();

    context_menu->addAction("Layer meta data", [&]() { showMetaData(true); });

    QMenu * settings_menu = new QMenu("Settings");
    settings_menu->addAction("Show/hide grid lines");
    settings_menu->addAction("Show/hide axis legends");
    context_menu->addSeparator();

    context_menu->addAction("Switch to 3D view", [&]() { emit showCurrentPeaksAs3D(); });

    const RangeType e_units = [&](){ // mouse position in units
      RangeType r;
      unit_mapper_.fromXY(widgetToData_(e->pos()), r);
      return r;
    }();

    // a small 10x10 pixel area around the current mouse position
    auto check_area = visible_area_.cloneWith({widgetToData_(e->pos() - QPoint(10, 10)), widgetToData_(e->pos() + QPoint(10, 10))}).getAreaUnit();

    //-------------------PEAKS----------------------------------
    if (auto* lp = dynamic_cast<const LayerDataPeak*>(&layer))
    {
      //add settings
      settings_menu->addSeparator();
      settings_menu->addAction("Show/hide projections");
      settings_menu->addAction("Show/hide MS/MS precursors");

      //add surrounding survey scans
      //find nearest survey scan
      SignedSize size = lp->getPeakData()->size();
      Int current = lp->getPeakData()->RTBegin(e_units.getMinRT()) - lp->getPeakData()->begin();
      if (current == size)  // if the user clicked right of the last MS1 scan
      {
        current = std::max(SignedSize{0}, size - 1); // we want the rightmost valid scan index
      }

      SignedSize i = 0;
      while (current + i < size || current - i >= 0)
      {
        if (current + i < size && (*lp->getPeakData())[current + i].getMSLevel() == 1)
        {
          current += i;
          break;
        }
        if (current - i >= 0 && (*lp->getPeakData())[current - i].getMSLevel() == 1)
        {
          current -= i;
          break;
        }
        ++i;
      }
      // search for four scans in both directions
      vector<Int> indices;
      indices.push_back(current);
      i = 1;
      while (current - i >= 0 && indices.size() < 5)
      {
        if ((*lp->getPeakData())[current - i].getMSLevel() == 1)
        {
          indices.push_back(current - i);
        }
        ++i;
      }
      i = 1;
      while (current + i < size && indices.size() < 9)
      {
        if ((*lp->getPeakData())[current + i].getMSLevel() == 1)
        {
          indices.push_back(current + i);
        }
        ++i;
      }
      sort(indices.rbegin(), indices.rend());
      QMenu* ms1_scans = context_menu->addMenu("Survey scan in 1D");
      QMenu* ms1_meta = context_menu->addMenu("Survey scan meta data");
      context_menu->addSeparator();
      for (i = 0; i < (Int)indices.size(); ++i)
      {
        if (indices[i] == current)
        {
          ms1_scans->addSeparator();
        }
        a = ms1_scans->addAction(QString("RT: ") + QString::number((*lp->getPeakData())[indices[i]].getRT()));
        a->setData(indices[i]);
        if (indices[i] == current)
        {
          ms1_scans->addSeparator();
        }

        if (indices[i] == current)
        {
          ms1_meta->addSeparator();
        }
        a = ms1_meta->addAction(QString("RT: ") + QString::number((*lp->getPeakData())[indices[i]].getRT()));
        a->setData(indices[i]);
        if (indices[i] == current)
        {
          ms1_meta->addSeparator();
        }
      }

      // add surrounding fragment scans
      // - We first attempt to look at the position where the user clicked
      // - Next we look within the +/- 5 scans around that position
      // - Next we look within the whole visible area
      QMenu* msn_scans = new QMenu("fragment scan in 1D");
      QMenu* msn_meta = new QMenu("fragment scan meta data");
      bool item_added = collectFragmentScansInArea_(check_area, a, msn_scans, msn_meta);
      if (!item_added)
      {
        // Now simply go for the 5 closest points in RT and check whether there
        // are any scans.
        // NOTE: that if we go for the visible area, we run the
        // risk of iterating through *all* the scans.
        check_area.RangeMZ::extend((RangeMZ)visible_area_.getAreaUnit());
        const auto& specs = lp->getPeakData()->getSpectra();
        check_area.RangeRT::operator=(RangeRT(specs[indices.back()].getRT(), specs[indices.front()].getRT()));
        item_added = collectFragmentScansInArea_(check_area, a, msn_scans, msn_meta);

        if (!item_added)
        { // OK, now lets search the whole visible area (may be large!)
          item_added = collectFragmentScansInArea_(visible_area_.getAreaUnit(), a, msn_scans, msn_meta);
        }
      }
      if (item_added)
      {
        context_menu->addMenu(msn_scans);
        context_menu->addMenu(msn_meta);
        context_menu->addSeparator();
      }
      
      auto it_closest_MS = lp->getPeakData()->getClosestSpectrumInRT(e_units.getMinRT());
      if (it_closest_MS->containsIMData())
      {
        context_menu->addAction(("Switch to ion mobility view (MSLevel: " + String(it_closest_MS->getMSLevel()) + ";RT: " + String(it_closest_MS->getRT(), false) + ")").c_str(),
                                [&]() {emit showCurrentPeaksAsIonMobility(*it_closest_MS); });
      }


      finishContextMenu_(context_menu, settings_menu);

      // evaluate menu
      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->parent() == ms1_scans  || result->parent() == msn_scans)
        {
          emit showSpectrumAsNew1D(result->data().toInt());
        }
        else if (result->parent() == ms1_meta || result->parent() == msn_meta)
        {
          showMetaData(true, result->data().toInt());
        }
      }
    }
    //-------------------FEATURES----------------------------------
    else if (auto* lf = dynamic_cast<const LayerDataFeature*>(&layer))
    {
      // add settings
      settings_menu->addSeparator();
      settings_menu->addAction("Show/hide convex hull");
      settings_menu->addAction("Show/hide trace convex hulls");
      settings_menu->addAction("Show/hide numbers/labels");
      settings_menu->addAction("Show/hide unassigned peptide hits");

      // search for nearby features
      QMenu* meta = new QMenu("Feature meta data");
      const FeatureMapType& features = *lf->getFeatureMap();
      // feature meta data menu
      for (auto it = features.cbegin(); it != features.cend(); ++it)
      {
        if (check_area.containsMZ(it->getMZ()) && check_area.containsRT(it->getRT()))
        {
          a = meta->addAction(QString("RT: ") + QString::number(it->getRT()) + "  m/z:" + QString::number(it->getMZ())
                              + "  charge:" + QString::number(it->getCharge()));
          a->setData((int)(it - features.begin()));
        }
      }
      if (! meta->actions().empty())
      {
        context_menu->addMenu(meta);
        context_menu->addSeparator();
      }

      // add modifiable flag
      settings_menu->addSeparator();
      settings_menu->addAction("Toggle edit/view mode", [&]() { getCurrentLayer().modifiable = ! getCurrentLayer().modifiable; });

      finishContextMenu_(context_menu, settings_menu);

      // evaluate menu
      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->text().left(3) == "RT:")
        {
          showMetaData(true, result->data().toInt());
        }
      }
    }
    //-------------------CONSENSUS FEATURES----------------------------------
    else if (auto* lc = dynamic_cast<const LayerDataConsensus*>(&layer))
    {
      // add settings
      settings_menu->addSeparator();
      settings_menu->addAction("Show/hide elements", [&]() { setLayerFlag(LayerDataBase::C_ELEMENTS, ! getLayerFlag(LayerDataBase::C_ELEMENTS)); });

      // search for nearby features
      QMenu* consens_meta = new QMenu("Consensus meta data");
      const ConsensusMapType& features = *lc->getConsensusMap();
      // consensus feature meta data menu
      for (auto it = features.cbegin(); it != features.cend(); ++it)
      {
        if (check_area.containsMZ(it->getMZ()) && check_area.containsRT(it->getRT()))
        {
          a = consens_meta->addAction(QString("RT: ") + QString::number(it->getRT()) + "  m/z:" + QString::number(it->getMZ()) + "  charge:" + QString::number(it->getCharge()));
          a->setData((int)(it - features.begin()));
        }
      }
      if (!consens_meta->actions().empty())
      {
        context_menu->addMenu(consens_meta);
        context_menu->addSeparator();
      }

      finishContextMenu_(context_menu, settings_menu);

      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->text().left(3) == "RT:")
        {
          showMetaData(true, result->data().toInt());
        }
      }
    }
    //------------------CHROMATOGRAMS----------------------------------
    else if (auto* lc = dynamic_cast<const LayerDataChrom*>(&layer))
    {
      settings_menu->addSeparator();
      settings_menu->addAction("Show/hide projections");
      settings_menu->addAction("Show/hide MS/MS precursors");

      const PeakMap& exp = *lc->getChromatogramData();

      constexpr int CHROMATOGRAM_SHOW_MZ_RANGE = 10;
      auto search_area = e_units;
      search_area.RangeMZ::extendLeftRight(CHROMATOGRAM_SHOW_MZ_RANGE);

      // collect all precursor that fall into the mz rt window
      typedef std::set<Precursor, Precursor::MZLess> PCSetType;
      PCSetType precursor_in_rt_mz_window;
      for (auto iter = exp.getChromatograms().cbegin(); iter != exp.getChromatograms().cend(); ++iter)
      {
        if (search_area.containsMZ(iter->getPrecursor().getMZ()) &&
            RangeBase{iter->front().getRT(), iter->back().getRT()}.contains(search_area.getRangeForDim(MSDim::MZ)))
        {
          precursor_in_rt_mz_window.insert(iter->getPrecursor());
        }
      }

      // determine product chromatograms for each precursor
      map<Precursor, vector<Size>, Precursor::MZLess> map_precursor_to_chrom_idx;
      for (PCSetType::const_iterator pit = precursor_in_rt_mz_window.begin(); pit != precursor_in_rt_mz_window.end(); ++pit)
      {
        for (vector<MSChromatogram >::const_iterator iter = exp.getChromatograms().begin(); iter != exp.getChromatograms().end(); ++iter)
        {
          if (iter->getPrecursor() == *pit)
          {
            map_precursor_to_chrom_idx[*pit].push_back(iter - exp.getChromatograms().begin());
          }
        }
      }

      QMenu* msn_chromatogram  = nullptr;
      QMenu* msn_chromatogram_meta = nullptr;

      if (!map_precursor_to_chrom_idx.empty())
      {
        msn_chromatogram = context_menu->addMenu("Chromatogram");
        msn_chromatogram_meta = context_menu->addMenu("Chromatogram meta data");
        context_menu->addSeparator();

        for (auto mit = map_precursor_to_chrom_idx.cbegin(); mit != map_precursor_to_chrom_idx.cend(); ++mit)
        {
          // Show the peptide sequence if available, otherwise show the m/z and charge only
          QString precursor_string = QString("Precursor m/z: (")  + String(mit->first.getCharge()).toQString() + ") " + QString::number(mit->first.getMZ());
          if (mit->first.metaValueExists("peptide_sequence"))
          {
            precursor_string = QString::number(mit->first.getMZ()) + " : " + String(mit->first.getMetaValue("peptide_sequence")).toQString() + " (" + QString::number(mit->first.getCharge()) + "+)";
          }
          QMenu * msn_precursor = msn_chromatogram->addMenu(precursor_string);  // new entry for every precursor

          // Show all: iterate over all chromatograms corresponding to the current precursor and add action containing all chromatograms
          a = msn_precursor->addAction(QString("Show all"));
          QList<QVariant> chroms_idx;
          for (auto vit = mit->second.cbegin(); vit != mit->second.cend(); ++vit)
          {
            chroms_idx.push_back((unsigned int)*vit);
          }
          a->setData(chroms_idx);

          // Show single chromatogram: iterate over all chromatograms corresponding to the current precursor and add action for the single chromatogram
          for (auto vit = mit->second.cbegin(); vit != mit->second.cend(); ++vit)
          {
            a = msn_precursor->addAction(QString("Chromatogram m/z: ") + QString::number(exp.getChromatograms()[*vit].getMZ())); // Precursor => Chromatogram MZ
            a->setData((int)(*vit));
          }
        }
      }

      finishContextMenu_(context_menu, settings_menu);

      // show context menu and evaluate result
      if ((result = context_menu->exec(mapToGlobal(e->pos()))))
      {
        if (result->parent()->parent() == msn_chromatogram) // clicked on chromatogram entry (level 2)
        {
          if (result->text() == "Show all")
          {
            std::vector<int> chrom_indices;
            for (const auto& var : result->data().toList())
            {
              chrom_indices.push_back(var.toInt());
              cout << "chrom_indices: " << var.toInt() << std::endl;
            }
            emit showChromatogramsAsNew1D(chrom_indices);
          }
          else   // Show single chromatogram
          {
            //cout << "Chromatogram result " << result->data().toInt() << endl;
            emit showSpectrumAsNew1D(result->data().toInt());
          }
        }
        else if (result->parent() == msn_chromatogram_meta)
        {
          showMetaData(true, result->data().toInt());
        }
      }
    }

    // common actions of peaks and features
    if (result)
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
      else if (result->text() == "Layer" || result->text() == "Visible layer data")
      {
        saveCurrentLayer(result->text() == "Visible layer data");
      }
      else if (result->text() == "As image")
      {
        spectrum_widget_->saveAsImage();
      }
      else if (result->text() == "Show/hide projections")
      {
        emit toggleProjections();
      }
      else if (result->text() == "Show/hide MS/MS precursors")
      {
        setLayerFlag(LayerDataBase::P_PRECURSORS, !getLayerFlag(LayerDataBase::P_PRECURSORS));
      }
      else if (result->text() == "Show/hide convex hull")
      {
        setLayerFlag(LayerDataBase::F_HULL, !getLayerFlag(LayerDataBase::F_HULL));
      }
      else if (result->text() == "Show/hide trace convex hulls")
      {
        setLayerFlag(LayerDataBase::F_HULLS, !getLayerFlag(LayerDataBase::F_HULLS));
      }
      else if (result->text() == "Show/hide unassigned peptide hits")
      {
        setLayerFlag(LayerDataBase::F_UNASSIGNED, !getLayerFlag(LayerDataBase::F_UNASSIGNED));
      }
      else if (result->text() == "Show/hide numbers/labels")
      {
        if (layer.label == LayerDataBase::L_NONE)
        {
          getCurrentLayer().label = LayerDataBase::L_META_LABEL;
        }
        else
        {
          getCurrentLayer().label = LayerDataBase::L_NONE;
        }
      }

    }

    e->accept();
  }

  void Plot2DCanvas::finishContextMenu_(QMenu* context_menu, QMenu* settings_menu)
  {
    // finish settings menu
    settings_menu->addSeparator();
    settings_menu->addAction("Preferences");

    // create save menu
    QMenu* save_menu = new QMenu("Save");
    save_menu->addAction("Layer");
    save_menu->addAction("Visible layer data");
    save_menu->addAction("As image");

    //add settings menu
    context_menu->addMenu(save_menu);
    context_menu->addMenu(settings_menu);

    //add external context menu
    if (context_add_)
    {
      context_menu->addSeparator();
      context_menu->addMenu(context_add_);
    }
  }

  void Plot2DCanvas::showCurrentLayerPreferences()
  {
    Internal::Plot2DPrefDialog dlg(this);
    LayerDataBase& layer = getCurrentLayer();

    ColorSelector * bg_color = dlg.findChild<ColorSelector *>("bg_color");
    MultiGradientSelector * gradient = dlg.findChild<MultiGradientSelector *>("gradient");
    QComboBox * feature_icon = dlg.findChild<QComboBox *>("feature_icon");
    QSpinBox * feature_icon_size = dlg.findChild<QSpinBox *>("feature_icon_size");

    bg_color->setColor(QColor(String(param_.getValue("background_color").toString()).toQString()));
    gradient->gradient().fromString(layer.param.getValue("dot:gradient"));
    feature_icon->setCurrentIndex(feature_icon->findText(String(layer.param.getValue("dot:feature_icon").toString()).toQString()));
    feature_icon_size->setValue((int)layer.param.getValue("dot:feature_icon_size"));

    if (dlg.exec())
    {
      param_.setValue("background_color", bg_color->getColor().name().toStdString());
      layer.param.setValue("dot:feature_icon", feature_icon->currentText().toStdString());
      layer.param.setValue("dot:feature_icon_size", feature_icon_size->value());
      layer.param.setValue("dot:gradient", gradient->gradient().toString());
      emit preferencesChange();
    }
  }

  void Plot2DCanvas::currentLayerParametersChanged_()
  {
    recalculateDotGradient_(getCurrentLayerIndex());

    update_buffer_ = true;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void Plot2DCanvas::updateLayer(Size i)
  {
    //update nearest peak
    selected_peak_.clear();
    recalculateRanges_();
    resetZoom(false);     //no repaint as this is done in intensityModeChange_() anyway
    intensityModeChange_();
    modificationStatus_(i, false);
  }

  void Plot2DCanvas::translateVisibleArea_(double x_axis_rel, double y_axis_rel)
  {
    auto xy = visible_area_.getAreaXY();
    const auto shift = xy.diagonal() * AreaXYType::PositionType{x_axis_rel, y_axis_rel};
    // publish (bounds checking is done inside changeVisibleArea_)
    changeVisibleArea_(visible_area_.cloneWith(xy + shift));
  }

  void Plot2DCanvas::translateLeft_(Qt::KeyboardModifiers /*m*/)
  {
    translateVisibleArea_( -0.05, 0.0 );
  }

  void Plot2DCanvas::translateRight_(Qt::KeyboardModifiers /*m*/)
  {
    translateVisibleArea_( 0.05, 0.0 );
  }

  void Plot2DCanvas::translateForward_()
  {
    translateVisibleArea_( 0.0, 0.05 );
  }

  void Plot2DCanvas::translateBackward_()
  {
    translateVisibleArea_(0.0, -0.05);
  }

  void Plot2DCanvas::keyPressEvent(QKeyEvent * e)
  {
    // CTRL+ALT (exactly, not e.g. CTRL+ALT+KEYPAD<X>|SHIFT...)
    // note that Qt::KeypadModifier is also a modifier which gets activated when any keypad key is pressed
    if (e->modifiers() == (Qt::ControlModifier | Qt::AltModifier))
    {
      String status_changed;
      // +Home (MacOSX small keyboard: Fn+ArrowLeft) => increase point size
      if ((e->key() == Qt::Key_Home) && (pen_size_max_ < PEN_SIZE_MAX_LIMIT))
      {
        ++pen_size_max_;
        status_changed = "Max. dot size increased to '" + String(pen_size_max_) + "'";
      }
      // +End (MacOSX small keyboard: Fn+ArrowRight) => decrease point size
      else if ((e->key() == Qt::Key_End) && (pen_size_max_ > PEN_SIZE_MIN_LIMIT))
      {
        --pen_size_max_;
        status_changed = "Max. dot size decreased to '" + String(pen_size_max_) + "'";
      }
      // +PageUp => increase min. coverage threshold
      else if (e->key() == Qt::Key_PageUp && canvas_coverage_min_ < CANVAS_COVERAGE_MIN_LIMITHIGH)
      {
        canvas_coverage_min_ += 0.05; // 5% steps
        status_changed = "Min. coverage threshold increased to '" + String(canvas_coverage_min_) + "'";
      }
      // +PageDown => decrease min. coverage threshold
      else if (e->key() == Qt::Key_PageDown && canvas_coverage_min_ > CANVAS_COVERAGE_MIN_LIMITLOW)
      {
        canvas_coverage_min_ -= 0.05; // 5% steps
        status_changed = "Min. coverage threshold decreased to '" + String(canvas_coverage_min_) + "'";
      }
      if (!status_changed.empty())
      {
        emit sendStatusMessage(status_changed, 0);
        update_buffer_ = true; // full repaint
        update_(OPENMS_PRETTY_FUNCTION); // schedule repaint
        return;
      }
    }

    // Delete features
    LayerDataBase& layer = getCurrentLayer();
    auto* lf = dynamic_cast<LayerDataFeature*>(&layer);
    if (e->key() == Qt::Key_Delete && getCurrentLayer().modifiable && lf && selected_peak_.isValid())
    {
      lf->getFeatureMap()->erase(lf->getFeatureMap()->begin() + selected_peak_.peak);
      selected_peak_.clear();
      update_buffer_ = true;
      update_(OPENMS_PRETTY_FUNCTION);
      modificationStatus_(getCurrentLayerIndex(), true);
      return;
    }

    // call parent class
    PlotCanvas::keyPressEvent(e);
  }

  void Plot2DCanvas::keyReleaseEvent(QKeyEvent * e)
  {
    //zoom if in zoom mode and a valid rectangle is selected
    if (action_mode_ == AM_ZOOM && rubber_band_.isVisible())
    {
      rubber_band_.hide();
      QRect rect = rubber_band_.geometry();
      if (rect.width() != 0 && rect.height() != 0)
      {
        AreaXYType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
        changeVisibleArea_(visible_area_.cloneWith(area), true, true);
      }
    }
    else if (action_mode_ == AM_MEASURE)
    {
      measurement_start_.clear();
      update_(OPENMS_PRETTY_FUNCTION);
    }

    // do the normal stuff
    PlotCanvas::keyReleaseEvent(e);
  }

  void Plot2DCanvas::mouseDoubleClickEvent(QMouseEvent * e)
  {
    LayerDataBase& current_layer = getCurrentLayer();
    auto* lf = dynamic_cast<LayerDataFeature*>(&current_layer);

    if (current_layer.modifiable && lf)
    {
      Feature tmp;
      if (selected_peak_.isValid())       //edit existing feature
      {
        FeatureEditDialog dialog(this);
        dialog.setFeature((*lf->getFeatureMap())[selected_peak_.peak]);
        if (dialog.exec())
        {
          tmp = dialog.getFeature();
          (*lf->getFeatureMap())[selected_peak_.peak] = tmp;
        }
      }
      else       //create new feature
      {
        tmp.setRT(widgetToData_(e->pos())[1]);
        tmp.setMZ(widgetToData_(e->pos())[0]);
        FeatureEditDialog dialog(this);
        dialog.setFeature(tmp);
        if (dialog.exec())
        {
          tmp = dialog.getFeature();
          lf->getFeatureMap()->push_back(tmp);
        }
      }

      // update gradient if the min/max intensity changes
      if (!lf->getFeatureMap()->getRange().containsIntensity(tmp.getIntensity()))
      {
        lf->getFeatureMap()->updateRanges();
        recalculateRanges_();
        intensityModeChange_();
      }
      else // just repaint to show the changes
      {
        update_buffer_ = true;
        update_(OPENMS_PRETTY_FUNCTION);
      }

      modificationStatus_(getCurrentLayerIndex(), true);
    }
  }

  void Plot2DCanvas::mergeIntoLayer(Size i, const FeatureMapSharedPtrType& map)
  {
    auto& layer = dynamic_cast<LayerDataFeature&>(layers_.getLayer(i));

    //reserve enough space
    layer.getFeatureMap()->reserve(layer.getFeatureMap()->size() + map->size());
    //add features
    for (Size j = 0; j < map->size(); ++j)
    {
      layer.getFeatureMap()->push_back((*map)[j]);
    }
    // update the layer and overall ranges (if necessary)
    auto old_range = layer.getFeatureMap()->getRange();
    layer.getFeatureMap()->updateRanges();
    if (!old_range.containsIntensity(layer.getFeatureMap()->getRangeForDim(MSDim::INT)))
    {
      intensityModeChange_();
    }
    // clear intensity range and compare the remaining dimensions
    old_range.RangeIntensity::clear();
    if (!old_range.containsAll(layer.getFeatureMap()->getRange()))
    {
      recalculateRanges_();
      resetZoom(true);
    }
  }

  void Plot2DCanvas::mergeIntoLayer(Size i, const ConsensusMapSharedPtrType& map)
  {
    auto& layer = dynamic_cast<LayerDataConsensus&>(layers_.getLayer(i));
    OPENMS_PRECONDITION(layer.type == LayerDataBase::DT_CONSENSUS, "Plot2DCanvas::mergeIntoLayer(i, map) non-consensus-feature layer selected");
    //reserve enough space
    layer.getConsensusMap()->reserve(layer.getConsensusMap()->size() + map->size());
    //add features
    for (Size j = 0; j < map->size(); ++j)
    {
      layer.getConsensusMap()->push_back((*map)[j]);
    }
    // update the layer and overall ranges (if necessary)
    auto old_range = layer.getConsensusMap()->getRange();
    layer.getConsensusMap()->updateRanges();
    if (!old_range.containsIntensity(layer.getConsensusMap()->getRangeForDim(MSDim::INT)))
    {
      intensityModeChange_();
    }
    // clear intensity range and compare the remaining dimensions
    old_range.RangeIntensity::clear();
    if (!old_range.containsAll(layer.getConsensusMap()->getRange()))
    {
      recalculateRanges_();
      resetZoom(true);
    }
  }

  void Plot2DCanvas::mergeIntoLayer(Size i, vector<PeptideIdentification> & peptides)
  {
    LayerDataBase& layer = layers_.getLayer(i);
    OPENMS_PRECONDITION(layer.type == LayerDataBase::DT_IDENT, "Plot2DCanvas::mergeIntoLayer(i, peptides) non-identification layer selected");
    
    auto& layer_peptides = dynamic_cast<IPeptideIds*>(&layer)->getPeptideIds();
    // reserve enough space
    layer_peptides.reserve(layer_peptides.size() + peptides.size());
    // insert peptides
    layer_peptides.insert(layer_peptides.end(), peptides.begin(),
                               peptides.end());
    // update the layer and overall ranges
    recalculateRanges_();
    resetZoom(true);
  }

    bool Plot2DCanvas::collectFragmentScansInArea_(const RangeType& range, QAction* a, QMenu* msn_scans, QMenu* msn_meta)
    {
      auto& layer = dynamic_cast<LayerDataPeak&>(getCurrentLayer());
      bool item_added = false;
      const auto last_RT = layer.getPeakData()->RTEnd(range.getMaxRT());
      for (ExperimentType::ConstIterator it = layer.getPeakData()->RTBegin(range.getMinRT());
                                         it != last_RT; ++it)
      {
        if (it->getPrecursors().empty()) continue;

        double mz = it->getPrecursors()[0].getMZ();
        if (it->getMSLevel() > 1 && range.containsMZ(mz))
        {
          a = msn_scans->addAction(QString("RT: ") + QString::number(it->getRT()) + " mz: " + QString::number(mz));
          a->setData((int)(it - layer.getPeakData()->begin()));
          a = msn_meta->addAction(QString("RT: ") + QString::number(it->getRT()) + " mz: " + QString::number(mz));
          a->setData((int)(it - layer.getPeakData()->begin()));
          item_added = true;
        }
      }
      return item_added;
    }

} //namespace OpenMS
