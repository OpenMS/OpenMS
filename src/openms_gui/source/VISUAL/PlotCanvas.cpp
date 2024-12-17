// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimator.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/LayerData1DChrom.h>
#include <OpenMS/VISUAL/LayerData1DPeak.h>
#include <OpenMS/VISUAL/LayerDataChrom.h>
#include <OpenMS/VISUAL/LayerDataConsensus.h>
#include <OpenMS/VISUAL/LayerDataFeature.h>
#include <OpenMS/VISUAL/LayerDataIdent.h>
#include <OpenMS/VISUAL/LayerDataPeak.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>
#include <OpenMS/VISUAL/PlotCanvas.h>
#include <OpenMS/VISUAL/PlotWidget.h>
#include <OpenMS/VISUAL/VISITORS/LayerStoreData.h>

// QT
#include <QPaintEvent>
#include <QPainter>
#include <QtWidgets/QMessageBox>
#include <iostream>
#include <utility>

using namespace std;

namespace OpenMS
{
  PlotCanvas::PlotCanvas(const Param& /*preferences*/, QWidget* parent)
    : QWidget(parent),
      DefaultParamHandler("PlotCanvas"),
      unit_mapper_({DIM_UNIT::RT, DIM_UNIT::MZ}),
      visible_area_(&unit_mapper_),
      rubber_band_(QRubberBand::Rectangle, this)
  {
    // Prevent filling background
    setAttribute(Qt::WA_OpaquePaintEvent);
    // get mouse coordinates while mouse moves over diagram and for focus handling
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);

    setMinimumSize(200, 200);
    setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);

    // set common defaults for all canvases
    defaults_.setValue("default_path", ".", "Default path for loading/storing data.");

    // Set focus policy in order to get keyboard events

    // Set 'whats this' text
    setWhatsThis(
      "Translate: Translate mode is activated by default. Hold down the left mouse key and move the mouse to translate. Arrow keys can be used for translation independent of the current mode.\n\n"
      "Zoom: Zoom mode is activated with the CTRL key. CTRL+/CTRL- are used to traverse the zoom stack (or mouse wheel). Pressing Backspace resets the zoom.\n\n"
      "Measure: Measure mode is activated with the SHIFT key. To measure the distance between data points, press the left mouse button on a point and drag the mouse to another point.\n\n");

    // set move cursor and connect signal that updates the cursor automatically
    updateCursor_();
    connect(this, SIGNAL(actionModeChange()), this, SLOT(updateCursor_()));
  }

  PlotCanvas::~PlotCanvas() = default;

  void PlotCanvas::resizeEvent(QResizeEvent* /* e */)
  {
#ifdef DEBUG_TOPPVIEW
    cout << "BEGIN " << OPENMS_PRETTY_FUNCTION << endl;
#endif
    buffer_ = QImage(width(), height(), QImage::Format_RGB32);
    update_buffer_ = true;
    updateScrollbars_();
    update_(OPENMS_PRETTY_FUNCTION);
#ifdef DEBUG_TOPPVIEW
    cout << "END   " << OPENMS_PRETTY_FUNCTION << endl;
#endif
  }

  void PlotCanvas::setFilters(const DataFilters& filters)
  {
    // set filters
    layers_.getCurrentLayer().filters = filters;
    // update the content
    update_buffer_ = true;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void PlotCanvas::showGridLines(bool show)
  {
    show_grid_ = show;
    update_buffer_ = true;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void PlotCanvas::intensityModeChange_()
  {
    // update axes (e.g. make it Log-scale)
    if (spectrum_widget_)
    {
      spectrum_widget_->updateAxes();
    }
    recalculateSnapFactor_();
    update_buffer_ = true;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void PlotCanvas::dimensionsChanged_()
  {
    zoom_stack_.clear(); // any zoom history is bogus

    // swap axes if necessary
    if (spectrum_widget_)
    {
      spectrum_widget_->updateAxes();
    }

    updateScrollbars_();
    update_buffer_ = true;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void PlotCanvas::changeVisibleArea_(VisibleArea new_area, bool repaint, bool add_to_stack)
  {
    auto data_range = getDataRange(); // getDataRange() is virtual, since its special for 1D (0-based intensity)
    if (intensity_mode_ == IM_PERCENTAGE)
    { // new_area will have [0, 100], and we don't want to make that any smaller if the data only goes up to, say 50
    }
    else
    { // make sure we stay inside the overall data range
      new_area.pushInto(data_range);
    }

    // store old zoom state
    if (add_to_stack)
    {
      // if we scrolled in between zooming we want to store the last position before zooming as well
      if ((!zoom_stack_.empty()) && (zoom_stack_.back() != visible_area_))
      {
        zoomAdd_(visible_area_);
      }
      // add current zoom
      zoomAdd_(new_area);
    }

    // always update, even if the area did not change, since the intensity mode might have changed
    visible_area_ = new_area;
    updateScrollbars_();
    recalculateSnapFactor_();
    emit visibleAreaChanged(new_area); // calls PlotWidget::updateAxes, which calls Plot(1D/2D/3D)Widget::recalculateAxes_
    emit layerZoomChanged(this); // calls TOPPViewBase::zoomOtherWindows (for linked windows)

    if (repaint)
    {
      update_buffer_ = true;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  void PlotCanvas::updateScrollbars_()
  {
  }

  void PlotCanvas::wheelEvent(QWheelEvent* e)
  {
    zoom_(e->position().x(), e->position().y(), e->angleDelta().y() > 0);
    e->accept();
  }

  void PlotCanvas::zoom_(int x, int y, bool zoom_in)
  {
    if (!zoom_in)
    {
      zoomBack_();
    }
    else
    { // we want to zoom into (x,y), which is in pixel units, hence we need to know the relative position of (x,y) in the widget
      constexpr PointXYType::CoordinateType zoom_factor = 0.8;
      const double rel_pos_x = (PointXYType::CoordinateType)x / width();
      const double rel_pos_y = (PointXYType::CoordinateType)(height() - y) / height();
      auto new_area = visible_area_.getAreaXY();
      {
        auto zoomed = Math::zoomIn(new_area.minX(), new_area.maxX(), zoom_factor, rel_pos_x);
        new_area.setMinX(zoomed.first);
        new_area.setMaxX(zoomed.second);
      }
      {
        auto zoomed = Math::zoomIn(new_area.minY(), new_area.maxY(), zoom_factor, rel_pos_y);
        new_area.setMinY(zoomed.first);
        new_area.setMaxY(zoomed.second);
      }

      if (new_area != visible_area_.getAreaXY())
      {
        zoomAdd_(visible_area_.cloneWith(new_area));
        PlotCanvas::changeVisibleArea_(*zoom_pos_);
      }
    }
  }

  void PlotCanvas::zoomBack_()
  {
    if (zoom_pos_ != zoom_stack_.begin())
    {
      --zoom_pos_;
      changeVisibleArea_(*zoom_pos_);
    }
  }

  void PlotCanvas::zoomForward_()
  {
    // if at end of zoom level then simply add a new zoom
    if (zoom_pos_ == zoom_stack_.end() || (zoom_pos_ + 1) == zoom_stack_.end())
    {
      auto new_area = visible_area_;
      auto xy = new_area.getAreaXY();
      zoomAdd_(new_area.setArea(xy.extend(0.8)));
      zoom_pos_ = --zoom_stack_.end(); // set to last position
    }
    else // goto next zoom level
    {
      ++zoom_pos_;
    }
    changeVisibleArea_(*zoom_pos_);
  }

  void PlotCanvas::zoomAdd_(const VisibleArea& area)
  {
    if (zoom_pos_ != zoom_stack_.end() && (zoom_pos_ + 1) != zoom_stack_.end())
    {
      zoom_stack_.erase(zoom_pos_ + 1, zoom_stack_.end());
    }
    zoom_stack_.push_back(area);
    zoom_pos_ = zoom_stack_.end();
    --zoom_pos_;
  }

  void PlotCanvas::zoomClear_()
  {
    zoom_stack_.clear();
    zoom_pos_ = zoom_stack_.end();
  }

  void PlotCanvas::resetZoom(bool repaint)
  {
    zoomClear_();
    changeVisibleArea_(visible_area_.cloneWith(overall_data_range_), repaint, true);
  }

  void PlotCanvas::setVisibleArea(const VisibleArea& area)
  { // do not simply call "changeVisibleArea_(area);", since this will choke on different
    // internal DimMappers (and you probably do not want to change the DimMapping. E.g. when calling this from a 2DCanvas (RT,mz) to display a 1DCanvas (mz,int))
    changeVisibleArea_(visible_area_.cloneWith(area.getAreaUnit()));
  }

  void PlotCanvas::setVisibleArea(const RangeAllType& area)
  {
    changeVisibleArea_(visible_area_.cloneWith(area));
  }

  void PlotCanvas::setVisibleArea(const AreaXYType& area)
  {
    changeVisibleArea_(visible_area_.cloneWith(area));
  }

  void PlotCanvas::setVisibleAreaX(double min, double max)
  {
    auto va = visible_area_.getAreaXY();
    va.setMinX(min);
    va.setMaxX(max);
    setVisibleArea(va);
  }

  void PlotCanvas::setVisibleAreaY(double min, double max)
  {
    auto va = visible_area_.getAreaXY();
    va.setMinY(min);
    va.setMaxY(max);
    setVisibleArea(va);
  }


  void PlotCanvas::saveCurrentLayer(bool visible)
  {
    const LayerDataBase& layer = getCurrentLayer();

    // determine proposed filename
    String proposed_name = param_.getValue("default_path").toString();
    if (visible == false && layer.filename != "")
    {
      proposed_name = layer.filename;
    }

    auto formats = layer.storeFullData()->getSupportedFileFormats(); // storeFullData() is cheap; we just want the formats...
    QString file_name = GUIHelpers::getSaveFilename(this, "Save file", proposed_name.toQString(), formats, true, formats.getTypes().front());
    if (file_name.isEmpty())
    {
      return;
    }

    auto visitor_data = visible ? layer.storeVisibleData(getVisibleArea().getAreaUnit(), layer.filters) : layer.storeFullData();
    visitor_data->saveToFile(file_name, ProgressLogger::GUI);
    modificationStatus_(getCurrentLayerIndex(), false);
  }

  void PlotCanvas::paintGridLines_(QPainter& painter)
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

    unsigned int xl, xh, yl, yh; // width/height of the diagram area, x, y coordinates of lo/hi x,y values

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
        case 0: // style settings for big intervals
          painter.setPen(p1);
          break;

        case 1: // style settings for small intervals
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
        case 0: // style settings for big intervals
          painter.setPen(p1);
          break;

        case 1: // style settings for small intervals
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

        painter.drawLine(xl, y, xh, y);
      }
    }

    painter.restore();
  }


  void setBaseLayerParameters(LayerDataBase* new_layer, const Param& param, const String& filename, const String& caption)
  {
    new_layer->param = param;
    new_layer->filename = filename;
    if (! caption.empty())
    {
      new_layer->setName(caption);
    }
    else 
    {
      new_layer->setName(QFileInfo(filename.toQString()).completeBaseName());
    }
  }

  bool PlotCanvas::addLayer(std::unique_ptr<LayerData1DBase> new_layer)
  {
    setBaseLayerParameters(new_layer.get(), param_, new_layer->filename, new_layer->getName());
    layers_.addLayer(std::move(new_layer));

    return finishAdding_();
  }

  bool PlotCanvas::addPeakLayer(const ExperimentSharedPtrType& map,
                                ODExperimentSharedPtrType od_map,
                                const String& filename,
                                const String& caption,
                                const bool use_noise_cutoff)
  {
    if (map->getSpectra().empty())
    {
      auto msg = "Your input data contains no spectra. Not adding layer.";
      OPENMS_LOG_WARN << msg << std::endl;
      QMessageBox::critical(this, "Error", msg);
      return false;
    }

    LayerDataPeakUPtr new_layer;
    if (dynamic_cast<Plot1DCanvas*>(this))
      new_layer.reset(new LayerData1DPeak);
    else
      new_layer.reset(new LayerDataPeak);
    new_layer->setPeakData(map);
    new_layer->setOnDiscPeakData(std::move(od_map));

    setBaseLayerParameters(new_layer.get(), param_, filename, caption);
    layers_.addLayer(std::move(new_layer));

    // calculate noise
    if (use_noise_cutoff)
    {
      auto cutoff = estimateNoiseFromRandomScans(*map, 1, 10, 5); // 5% of low intensity data is considered noise
      DataFilters filters;
      filters.add(DataFilters::DataFilter(DataFilters::INTENSITY, DataFilters::GREATER_EQUAL, cutoff));
      setFilters(filters);
    }
    else // no mower, hide zeros if wanted
    {
      if (map->hasZeroIntensities(1))
      {
        DataFilters filters;
        filters.add(DataFilters::DataFilter(DataFilters::INTENSITY, DataFilters::GREATER_EQUAL, 0.001));
        setFilters(filters);
      }
    }

    return finishAdding_();
  }

  
  bool PlotCanvas::addChromLayer(const ExperimentSharedPtrType& map, ODExperimentSharedPtrType od_map, const String& filename, const String& caption)
  {
    if (map->getChromatograms().empty())
    {
      auto msg = "Your input data contains no chromatograms. Not adding layer.";
      OPENMS_LOG_WARN << msg << std::endl;
      QMessageBox::critical(this, "Error", msg);
      return false;
    }

    LayerDataChromUPtr new_layer;
    if (dynamic_cast<Plot1DCanvas*>(this))
      new_layer.reset(new LayerData1DChrom);
    else
      new_layer.reset(new LayerDataChrom);
    new_layer->setChromData(map);
    new_layer->setOnDiscPeakData(std::move(od_map));

    setBaseLayerParameters(new_layer.get(), param_, filename, caption);
    layers_.addLayer(std::move(new_layer));

    return finishAdding_();
  }

  bool PlotCanvas::addLayer(FeatureMapSharedPtrType map, const String& filename, const String& caption)
  {
    LayerDataFeatureUPtr new_layer(new LayerDataFeature);
    new_layer->getFeatureMap() = std::move(map);

    setBaseLayerParameters(new_layer.get(), param_, filename, caption);
    layers_.addLayer(std::move(new_layer));
    return finishAdding_();
  }

  bool PlotCanvas::addLayer(ConsensusMapSharedPtrType map, const String& filename, const String& caption)
  {
    LayerDataBaseUPtr new_layer(new LayerDataConsensus(map));

    setBaseLayerParameters(new_layer.get(), param_, filename, caption);
    layers_.addLayer(std::move(new_layer));
    return finishAdding_();
  }

  bool PlotCanvas::addLayer(vector<PeptideIdentification>& peptides, const String& filename, const String& caption)
  {
    LayerDataIdent* new_layer(new LayerDataIdent); // ownership will be transferred to unique_ptr below; no need to delete
    new_layer->setPeptideIds(std::move(peptides));
    setBaseLayerParameters(new_layer, param_, filename, caption);
    layers_.addLayer(LayerDataBaseUPtr(new_layer));
    return finishAdding_();
  }

  void PlotCanvas::popIncompleteLayer_(const QString& error_message)
  {
    layers_.removeCurrentLayer();
    if (!error_message.isEmpty())
      QMessageBox::critical(this, "Error", error_message);
  }

  void PlotCanvas::setLayerName(Size i, const String& name)
  {
    getLayer(i).setName(name);
    if (i == 0 && spectrum_widget_)
    {
      spectrum_widget_->setWindowTitle(name.toQString());
    }
  }

  String PlotCanvas::getLayerName(const Size i)
  {
    return getLayer(i).getName();
  }

  void PlotCanvas::changeVisibility(Size i, bool b)
  {
    LayerDataBase& layer = getLayer(i);
    if (layer.visible != b)
    {
      layer.visible = b;
      update_buffer_ = true;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  void PlotCanvas::changeLayerFilterState(Size i, bool b)
  {
    LayerDataBase& layer = getLayer(i);
    if (layer.filters.isActive() != b)
    {
      layer.filters.setActive(b);
      update_buffer_ = true;
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  const PlotCanvas::RangeType& PlotCanvas::getDataRange() const
  {
    return overall_data_range_;
  }

  void PlotCanvas::recalculateRanges_()
  {
    RangeType& layer_range = overall_data_range_;
    layer_range.clearRanges();

    for (Size layer_index = 0; layer_index < getLayerCount(); ++layer_index)
    {
      layer_range.extend(getLayer(layer_index).getRange());
    }
    // set minimum intensity to 0 (avoid negative intensities!)
    if (layer_range.getMinIntensity() < 0) layer_range.setMinIntensity(0);

    // add 4% margin (2% left, 2% right) to RT, m/z, IM and intensity
    layer_range.scaleBy(1.04);

    // make sure that each dimension is not a single point (axis widget won't like that)
    // (this needs to be the last command to ensure this property holds when leaving the function!)
    layer_range.minSpanIfSingular(1);
  }

  double PlotCanvas::getSnapFactor()
  { // only useful for 1D view at the moment (which only has a single snap factor). 2D view has as many as there are layers.
    return snap_factors_[0];
  }

  double PlotCanvas::getPercentageFactor() const
  {
    return percentage_factor_;
  }

  void PlotCanvas::recalculateSnapFactor_()
  {
  }

  void PlotCanvas::horizontalScrollBarChange(int /*value*/)
  {
  }

  void PlotCanvas::verticalScrollBarChange(int /*value*/)
  {
  }

  void PlotCanvas::update_(const char*)
  {
    update();
  }

  // this does not work anymore, probably due to Qt::StrongFocus :( -- todo: delete!
  void PlotCanvas::focusOutEvent(QFocusEvent* /*e*/)
  {
    // Alt/Shift pressed and focus lost => change back action mode
    if (action_mode_ != AM_TRANSLATE)
    {
      action_mode_ = AM_TRANSLATE;
      emit actionModeChange();
    }

    // reset peaks
    selected_peak_.clear();
    measurement_start_.clear();

    // update
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void PlotCanvas::leaveEvent(QEvent* /*e*/)
  {
    // release keyboard, when the mouse pointer leaves
    releaseKeyboard();
  }

  void PlotCanvas::enterEvent(QEnterEvent* /*e*/)
  {
    // grab keyboard, as we need to handle key presses
    grabKeyboard();
  }

  void PlotCanvas::keyReleaseEvent(QKeyEvent* e)
  {
    // Alt/Shift released => change back action mode
    if (e->key() == Qt::Key_Control || e->key() == Qt::Key_Shift)
    {
      action_mode_ = AM_TRANSLATE;
      emit actionModeChange();
      e->accept();
    }

    e->ignore();
  }

  void PlotCanvas::keyPressEvent(QKeyEvent* e)
  {
    if (e->key() == Qt::Key_Control)
    { // Ctrl pressed => change action mode
      action_mode_ = AM_ZOOM;
      emit actionModeChange();
    }
    else if (e->key() == Qt::Key_Shift)
    { // Shift pressed => change action mode
      action_mode_ = AM_MEASURE;
      emit actionModeChange();
    }
    else if ((e->modifiers() & Qt::ControlModifier) &&
             (e->key() == Qt::Key_Plus)) // do not use (e->modifiers() == Qt::ControlModifier) to target Ctrl exclusively, since +/- might(!) also trigger the Qt::KeypadModifier
    {                                    // CTRL+Plus => Zoom stack
      zoomForward_();
    }
    else if ((e->modifiers() & Qt::ControlModifier) && (e->key() == Qt::Key_Minus))
    { // CTRL+Minus => Zoom stack
      zoomBack_();
    }
    // Arrow keys => translate
    else if (e->key() == Qt::Key_Left)
    {
      translateLeft_(e->modifiers());
    }
    else if (e->key() == Qt::Key_Right)
    {
      translateRight_(e->modifiers());
    }
    else if (e->key() == Qt::Key_Up)
    {
      translateForward_();
    }
    else if (e->key() == Qt::Key_Down)
    {
      translateBackward_();
    }
    else if (e->key() == Qt::Key_Backspace)
    { // Backspace to reset zoom
      resetZoom();
    }
    else if ((e->modifiers() == (Qt::ControlModifier | Qt::AltModifier)) && (e->key() == Qt::Key_T))
    { // CTRL+ALT+T => activate timing mode
      show_timing_ = !show_timing_;
    }
    else
    { // call the keyPressEvent() of the parent widget
      e->ignore();
    }
  }

  void PlotCanvas::translateLeft_(Qt::KeyboardModifiers /*m*/)
  {
  }

  void PlotCanvas::translateRight_(Qt::KeyboardModifiers /*m*/)
  {
  }

  void PlotCanvas::translateForward_()
  {
  }

  void PlotCanvas::translateBackward_()
  {
  }

  void PlotCanvas::setAdditionalContextMenu(QMenu* menu)
  {
    context_add_ = menu;
  }

  const DimMapper<2>& PlotCanvas::getMapper() const
  {
    return unit_mapper_;
  }
  
  void PlotCanvas::setMapper(const DimMapper<2>& mapper)
  {
    unit_mapper_ = mapper;
  }

  void PlotCanvas::showMetaData(bool modifiable, Int index)
  {
    LayerDataBase& layer = getCurrentLayer();

    MetaDataBrowser dlg(modifiable, this);
    if (index == -1)
    {
      if (auto lp = dynamic_cast<LayerDataPeak*>(&layer))
      {
        dlg.add(*lp->getPeakDataMuteable());
        // Exception for Plot1DCanvas, here we add the meta data of the one spectrum
        if (auto lp1 = dynamic_cast<LayerData1DPeak*>(&layer))
        {
          dlg.add((*lp1->getPeakDataMuteable())[lp1->getCurrentIndex()]);
        }
      }
      if (auto lp = dynamic_cast<LayerDataFeature*>(&layer))
      {
        dlg.add(*lp->getFeatureMap());
      }
      if (auto lp = dynamic_cast<LayerDataConsensus*>(&layer))
      {
        dlg.add(*lp->getConsensusMap());
      }
      else if (layer.type == LayerDataBase::DT_CHROMATOGRAM)
      {
        // TODO CHROM
      }
      else if (layer.type == LayerDataBase::DT_IDENT)
      {
        // TODO IDENT
      }
    }
    else // show element meta data
    {
      if (auto lp = dynamic_cast<LayerDataPeak*>(&layer))
      {
        dlg.add((*lp->getPeakDataMuteable())[index]);
      }
      else if (auto lp = dynamic_cast<LayerDataFeature*>(&layer))
      {
        dlg.add((*lp->getFeatureMap())[index]);
      }
      else if (auto lp = dynamic_cast<LayerDataConsensus*>(&layer))
      {
        dlg.add((*lp->getConsensusMap())[index]);
      }
      else if (layer.type == LayerDataBase::DT_CHROMATOGRAM)
      {
        // TODO CHROM
      }
      else if (layer.type == LayerDataBase::DT_IDENT)
      {
        // TODO IDENT
      }
    }

    // if the meta data was modified, set the flag
    if (modifiable && dlg.exec())
    {
      modificationStatus_(getCurrentLayerIndex(), true);
    }
  }

  void PlotCanvas::updateCursor_()
  {
    switch (action_mode_)
    {
      case AM_TRANSLATE:
        setCursor(QCursor(QPixmap(":/cursor_move.png"), 0, 0));
        break;

      case AM_ZOOM:
        setCursor(QCursor(QPixmap(":/cursor_zoom.png"), 0, 0));
        break;

      case AM_MEASURE:
        setCursor(QCursor(QPixmap(":/cursor_measure.png"), 0, 0));
        break;
    }
  }

  void PlotCanvas::modificationStatus_(Size layer_index, bool modified)
  {
    LayerDataBase& layer = getLayer(layer_index);
    if (layer.modified != modified)
    {
      layer.modified = modified;
#ifdef DEBUG_TOPPVIEW
      cout << "BEGIN " << OPENMS_PRETTY_FUNCTION << endl;
      cout << "emit: layerModificationChange" << endl;
      cout << "END " << OPENMS_PRETTY_FUNCTION << endl;
#endif
      emit layerModficationChange(getCurrentLayerIndex(), modified);
    }
  }

  void PlotCanvas::drawText_(QPainter& painter, const QStringList& text)
  {
    GUIHelpers::drawText(painter, text, {2, 3}, Qt::black, QColor(255, 255, 255, 200));
  }

  double PlotCanvas::getIdentificationMZ_(const Size layer_index, const PeptideIdentification& peptide) const
  {
    if (getLayerFlag(layer_index, LayerDataBase::I_PEPTIDEMZ))
    {
      const PeptideHit& hit = peptide.getHits().front();
      Int charge = hit.getCharge();
      return hit.getSequence().getMZ(charge);
    }
    else
    {
      return peptide.getMZ();
    }
  }

  void LayerStack::addLayer(LayerDataBaseUPtr new_layer)
  {
    // insert after last layer of same type,
    // if there is no such layer after last layer of previous types,
    // if there are no layers at all put at front
    auto it = std::find_if(layers_.rbegin(), layers_.rend(), [&new_layer](const LayerDataBaseUPtr& l) { return l->type <= new_layer->type; });

    auto where = layers_.insert(it.base(), std::move(new_layer));
    // update to index we just inserted into
    current_layer_ = where - layers_.begin();
  }

  const LayerDataBase& LayerStack::getLayer(const Size index) const
  {
    if (index >= layers_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, layers_.size());
    }
    return *layers_[index].get();
  }

  LayerDataBase& LayerStack::getLayer(const Size index)
  {
    if (index >= layers_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, layers_.size());
    }
    return *layers_[index].get();
  }

  const LayerDataBase& LayerStack::getCurrentLayer() const
  {
    if (current_layer_ >= layers_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, current_layer_, layers_.size());
    }
    return *layers_[current_layer_].get();
  }

  LayerDataBase& LayerStack::getCurrentLayer()
  {
    if (current_layer_ >= layers_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, current_layer_, layers_.size());
    }
    return *layers_[current_layer_].get();
  }

  void LayerStack::setCurrentLayer(Size index)
  {
    if (index >= layers_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, layers_.size());
    }
    current_layer_ = index;
  }

  Size LayerStack::getCurrentLayerIndex() const
  {
    return current_layer_;
  }

  bool LayerStack::empty() const
  {
    return layers_.empty();
  }

  Size LayerStack::getLayerCount() const
  {
    return layers_.size();
  }

  void LayerStack::removeLayer(Size layer_index)
  {
    if (layer_index >= layers_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, layer_index, layers_.size());
    }
    layers_.erase(layers_.begin() + layer_index);

    // update current layer if it became invalid TODO dont you have to adjust the index to stay on the same layer??
    if (current_layer_ >= getLayerCount())
    {
      current_layer_ = getLayerCount() - 1; // overflow is intentional
    }
  }

  void LayerStack::removeCurrentLayer()
  {
    removeLayer(current_layer_);
  }

} // namespace OpenMS
