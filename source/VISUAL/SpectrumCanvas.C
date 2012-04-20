// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>

// QT
#include <QtGui/QPainter>
#include <QtGui/QPaintEvent>
#include <QtGui/QBitmap>
#include <QtGui/QWheelEvent>
#include <QtGui/QMessageBox>
#include <QtGui/QPushButton>
#include <QtGui/QFontMetrics>
#include <QtGui/QFontMetrics>

#include <iostream>

using namespace std;

namespace OpenMS
{
  SpectrumCanvas::SpectrumCanvas(const Param & /*preferences*/, QWidget * parent) :
    QWidget(parent),
    DefaultParamHandler("SpectrumCanvas"),
    buffer_(),
    action_mode_(AM_TRANSLATE),
    intensity_mode_(IM_NONE),
    layers_(),
    mz_to_x_axis_(true),
    visible_area_(AreaType::empty),
    overall_data_range_(DRange<3>::empty),
    show_grid_(true),
    zoom_stack_(),
    zoom_pos_(zoom_stack_.end()),
    update_buffer_(false),
    current_layer_(0),
    spectrum_widget_(0),
    percentage_factor_(1.0),
    snap_factors_(1, 1.0),
    rubber_band_(QRubberBand::Rectangle, this),
    context_add_(0),
    show_timing_(false),
    selected_peak_(),
    measurement_start_()
  {
    //Prevent filling background
    setAttribute(Qt::WA_OpaquePaintEvent);
    // get mouse coordinates while mouse moves over diagramm and for focus handling
    setMouseTracking(TRUE);
    setFocusPolicy(Qt::StrongFocus);

    setMinimumSize(200, 200);
    setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);

    //reserve enough space to avoid copying layer data
    layers_.reserve(10);

    //set common defaults for all canvases
    defaults_.setValue("default_path", ".", "Default path for loading/storing data.");

    //Set focus policy in order to get keyboard events

    //Set 'whats this' text
    setWhatsThis("Translate: Translate mode is activated by default. Hold down the left mouse key and move the mouse to translate. Arrow keys can be used for translation independent of the current mode.\n\n"
                 "Zoom: Zoom mode is activated with the CTRL key. CTRL+/CTRL- are used to traverse the zoom stack (or mouse wheel). Pressing Backspace resets the zoom.\n\n"
                 "Measure: Measure mode is activated with the SHIFT key. To measure the distace between data points, press the left mouse button on a point and drag the mouse to another point.\n\n"
                 );

    //set move cursor and connect signal that updates the cursor automatically
    updateCursor_();
    connect(this, SIGNAL(actionModeChange()), this, SLOT(updateCursor_()));
  }

  SpectrumCanvas::~SpectrumCanvas()
  {
    //cout << "DEST SpectrumCanvas" << endl;
  }

  void SpectrumCanvas::resizeEvent(QResizeEvent * /* e */)
  {
    #ifdef DEBUG_TOPPVIEW
    cout << "BEGIN " << __PRETTY_FUNCTION__ << endl;
    #endif
    buffer_ = QImage(width(), height(), QImage::Format_RGB32);
    update_buffer_ = true;
    updateScrollbars_();
    update_(__PRETTY_FUNCTION__);
    #ifdef DEBUG_TOPPVIEW
    cout << "END   " << __PRETTY_FUNCTION__ << endl;
    #endif
  }

  void SpectrumCanvas::setFilters(const DataFilters & filters)
  {
    //set filters
    layers_[current_layer_].filters = filters;
    //update the content
    update_buffer_ = true;
    update_(__PRETTY_FUNCTION__);
  }

  void SpectrumCanvas::showGridLines(bool show)
  {
    show_grid_ = show;
    update_buffer_ = true;
    update_(__PRETTY_FUNCTION__);
  }

  void SpectrumCanvas::intensityModeChange_()
  {
    recalculateSnapFactor_();
    update_buffer_ = true;
    update_(__PRETTY_FUNCTION__);
  }

  void SpectrumCanvas::mzToXAxis(bool mz_to_x_axis)
  {
    mz_to_x_axis_ = mz_to_x_axis;

    //swap axes if necessary
    if (spectrum_widget_)
    {
      spectrum_widget_->updateAxes();
    }

    updateScrollbars_();
    update_buffer_ = true;
    update_(__PRETTY_FUNCTION__);
  }

  void SpectrumCanvas::changeVisibleArea_(const AreaType & new_area, bool repaint, bool add_to_stack)
  {
    //store old zoom state
    if (add_to_stack)
    {
      // if we scrolled in between zooming we want to store the last position before zooming as well
      if ((zoom_stack_.size() > 0)
         && (zoom_stack_.back() != visible_area_))
      {
        zoomAdd_(visible_area_);
      }
      // add current zoom
      zoomAdd_(new_area);
    }

    if (new_area != visible_area_)
    {
      visible_area_ = new_area;
      updateScrollbars_();
      emit visibleAreaChanged(new_area);
    }

    if (repaint)
    {
      update_buffer_ = true;
      update_(__PRETTY_FUNCTION__);
    }
  }

  void SpectrumCanvas::updateScrollbars_()
  {

  }

  void SpectrumCanvas::wheelEvent(QWheelEvent * e)
  {
    if (e->delta() > 0)
    {
      zoomForward_();
    }
    else
    {
      zoomBack_();
    }
    e->accept();
  }

  void SpectrumCanvas::zoomBack_()
  {
    //cout << "Zoom out" << endl;
    //cout << " - pos before:" << (zoom_pos_-zoom_stack_.begin()) << endl;
    //cout << " - size before:" << zoom_stack_.size() << endl;
    if (zoom_pos_ != zoom_stack_.begin())
    {
      --zoom_pos_;
      changeVisibleArea_(*zoom_pos_);
    }
    //cout << " - pos after:" << (zoom_pos_-zoom_stack_.begin()) << endl;
  }

  void SpectrumCanvas::zoomForward_()
  {
    //cout << "Zoom in" << endl;
    //cout << " - pos before:" << (zoom_pos_-zoom_stack_.begin()) << endl;
    //cout << " - size before:" << zoom_stack_.size() <<endl;

    // if at end of zoom level then simply add a new zoom
    if (zoom_pos_ == zoom_stack_.end() || (zoom_pos_ + 1) == zoom_stack_.end())
    {
      AreaType new_area;
      // distance of areas center to border times a zoom factor of 0.8
      AreaType::CoordinateType size0 = visible_area_.width() / 2 * 0.8;
      AreaType::CoordinateType size1 = visible_area_.height() / 2 * 0.8;
      new_area.setMinX(visible_area_.center()[0] - size0);
      new_area.setMinY(visible_area_.center()[1] - size1);
      new_area.setMaxX(visible_area_.center()[0] + size0);
      new_area.setMaxY(visible_area_.center()[1] + size1);
      zoomAdd_(new_area);
      zoom_pos_ = --zoom_stack_.end(); // set to last position
    }
    else // goto next zoom level
    {
      ++zoom_pos_;
    }
    changeVisibleArea_(*zoom_pos_);

    //cout << " - pos after:" << (zoom_pos_-zoom_stack_.begin()) << endl;
  }

  void SpectrumCanvas::zoomAdd_(const AreaType & area)
  {
    //cout << "Adding to stack" << endl;
    //cout << " - pos before:" << (zoom_pos_-zoom_stack_.begin()) << endl;
    //cout << " - size before:" << zoom_stack_.size() <<endl;
    if (zoom_pos_ != zoom_stack_.end() && (zoom_pos_ + 1) != zoom_stack_.end())
    {
      //cout << " - removing from:" << ((zoom_pos_+1)-zoom_stack_.begin()) << endl;
      zoom_stack_.erase(zoom_pos_ + 1, zoom_stack_.end());
    }
    zoom_stack_.push_back(area);
    zoom_pos_ = zoom_stack_.end();
    --zoom_pos_;
    //cout << " - pos after:" << (zoom_pos_-zoom_stack_.begin()) << endl;
    //cout << " - size after:" << zoom_stack_.size() <<endl;
  }

  void SpectrumCanvas::zoomClear_()
  {
    zoom_stack_.clear();
    zoom_pos_ = zoom_stack_.end();
  }

  void SpectrumCanvas::resetZoom(bool repaint)
  {
    AreaType tmp;
    tmp.assign(overall_data_range_);
    zoomClear_();
    changeVisibleArea_(tmp, repaint, true);
  }

  void SpectrumCanvas::setVisibleArea(AreaType area)
  {
    //cout << __PRETTY_FUNCTION__ << endl;
    changeVisibleArea_(area);
  }

  void SpectrumCanvas::paintGridLines_(QPainter & painter)
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

      int x;
      for (std::vector<double>::const_iterator it = spectrum_widget_->xAxis()->gridLines()[j].begin(); it != spectrum_widget_->xAxis()->gridLines()[j].end(); ++it)
      {
        x = static_cast<int>(Math::intervalTransformation(*it, spectrum_widget_->xAxis()->getAxisMinimum(), spectrum_widget_->xAxis()->getAxisMaximum(), xl, xh));
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

      int y;
      for (std::vector<double>::const_iterator it = spectrum_widget_->yAxis()->gridLines()[j].begin(); it != spectrum_widget_->yAxis()->gridLines()[j].end(); ++it)
      {
        y = static_cast<int>(Math::intervalTransformation(*it, spectrum_widget_->yAxis()->getAxisMinimum(), spectrum_widget_->yAxis()->getAxisMaximum(), yl, yh));

        painter.drawLine(xl, y, xh, y);
      }
    }

    painter.restore();
  }

  Size SpectrumCanvas::activeLayerIndex() const
  {
    return current_layer_;
  }

  bool SpectrumCanvas::addLayer(ExperimentSharedPtrType map, const String & filename)
  {
    layers_.resize(layers_.size() + 1);
    layers_.back().param = param_;
    layers_.back().filename = filename;
    layers_.back().getPeakData() = map;

    if (layers_.back().getPeakData()->getChromatograms().size() != 0)
    {
      Size num_chrom(0);
      for (Size i = 0; i != layers_.back().getPeakData()->getChromatograms().size(); ++i)
      {
        if (layers_.back().getPeakData()->getChromatograms()[i].getChromatogramType() == ChromatogramSettings::SELECTED_ION_CURRENT_CHROMATOGRAM ||
            layers_.back().getPeakData()->getChromatograms()[i].getChromatogramType() == ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
        {
          ++num_chrom;
        }
      }

      if (num_chrom > 0)
      {
        layers_.back().type = LayerData::DT_CHROMATOGRAM;
      }
      else
      {
        layers_.back().type = LayerData::DT_PEAK;
      }
    }
    else
    {
      layers_.back().type = LayerData::DT_PEAK;
    }
    return finishAdding_();
  }

  bool SpectrumCanvas::addLayer(FeatureMapSharedPtrType map, const String & filename)
  {
    layers_.resize(layers_.size() + 1);
    layers_.back().param = param_;
    layers_.back().filename = filename;
    layers_.back().getFeatureMap() = map;
    layers_.back().type = LayerData::DT_FEATURE;
    return finishAdding_();
  }

  bool SpectrumCanvas::addLayer(ConsensusMapSharedPtrType map, const String & filename)
  {
    layers_.resize(layers_.size() + 1);
    layers_.back().param = param_;
    layers_.back().filename = filename;
    layers_.back().getConsensusMap() = map;
    layers_.back().type = LayerData::DT_CONSENSUS;
    return finishAdding_();
  }

  bool SpectrumCanvas::addLayer(vector<PeptideIdentification> & peptides,
                                const String & filename)
  {
    layers_.resize(layers_.size() + 1);
    layers_.back().param = param_;
    layers_.back().filename = filename;
    layers_.back().peptides.swap(peptides);
    layers_.back().type = LayerData::DT_IDENT;
    return finishAdding_();
  }

  void SpectrumCanvas::setLayerName(Size i, const String & name)
  {
    OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::setLayerName(i,name) index overflow");
    getLayer_(i).name = name;
    if (i == 0 && spectrum_widget_)
      spectrum_widget_->setWindowTitle(name.toQString());
  }

  String SpectrumCanvas::getLayerName(const Size i)
  {
    OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::getLayerName(i) index overflow");
    return getLayer_(i).name;
  }

  void SpectrumCanvas::changeVisibility(Size i, bool b)
  {
    OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::changeVisibility(i,b) index overflow");
    LayerData & layer = getLayer_(i);
    if (layer.visible != b)
    {
      layer.visible = b;
      update_buffer_ = true;
      update_(__PRETTY_FUNCTION__);
    }
  }

  void SpectrumCanvas::changeLayerFilterState(Size i, bool b)
  {
    OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::changeVisibility(i,b) index overflow");
    LayerData & layer = getLayer_(i);
    if (layer.filters.isActive() != b)
    {
      layer.filters.setActive(b);
      update_buffer_ = true;
      update_(__PRETTY_FUNCTION__);
    }
  }

  const DRange<3> & SpectrumCanvas::getDataRange()
  {
    return overall_data_range_;
  }

  void SpectrumCanvas::recalculateRanges_(UInt mz_dim, UInt rt_dim, UInt it_dim)
  {
    overall_data_range_ = DRange<3>::empty;
    DRange<3>::PositionType m_min = overall_data_range_.minPosition();
    DRange<3>::PositionType m_max = overall_data_range_.maxPosition();

    for (Size layer_index = 0; layer_index < getLayerCount(); ++layer_index)
    {
      if (getLayer(layer_index).type == LayerData::DT_PEAK || getLayer(layer_index).type == LayerData::DT_CHROMATOGRAM)
      {
        const ExperimentType & map = *getLayer(layer_index).getPeakData();
        if (map.getMinMZ() < m_min[mz_dim])
          m_min[mz_dim] = map.getMinMZ();
        if (map.getMaxMZ() > m_max[mz_dim])
          m_max[mz_dim] = map.getMaxMZ();
        if (map.getMinRT() < m_min[rt_dim])
          m_min[rt_dim] = map.getMinRT();
        if (map.getMaxRT() > m_max[rt_dim])
          m_max[rt_dim] = map.getMaxRT();
        if (map.getMinInt() < m_min[it_dim])
          m_min[it_dim] = map.getMinInt();
        if (map.getMaxInt() > m_max[it_dim])
          m_max[it_dim] = map.getMaxInt();
      }
      else if (getLayer(layer_index).type == LayerData::DT_FEATURE)
      {
        const FeatureMapType & map = *getLayer(layer_index).getFeatureMap();
        if (map.getMin()[1] < m_min[mz_dim])
          m_min[mz_dim] = map.getMin()[1];
        if (map.getMax()[1] > m_max[mz_dim])
          m_max[mz_dim] = map.getMax()[1];
        if (map.getMin()[0] < m_min[rt_dim])
          m_min[rt_dim] = map.getMin()[0];
        if (map.getMax()[0] > m_max[rt_dim])
          m_max[rt_dim] = map.getMax()[0];
        if (map.getMinInt() < m_min[it_dim])
          m_min[it_dim] = map.getMinInt();
        if (map.getMaxInt() > m_max[it_dim])
          m_max[it_dim] = map.getMaxInt();
      }
      else if (getLayer(layer_index).type == LayerData::DT_CONSENSUS)
      {
        const ConsensusMapType & map = *getLayer(layer_index).getConsensusMap();
        if (map.getMin()[1] < m_min[mz_dim])
          m_min[mz_dim] = map.getMin()[1];
        if (map.getMax()[1] > m_max[mz_dim])
          m_max[mz_dim] = map.getMax()[1];
        if (map.getMin()[0] < m_min[rt_dim])
          m_min[rt_dim] = map.getMin()[0];
        if (map.getMax()[0] > m_max[rt_dim])
          m_max[rt_dim] = map.getMax()[0];
        if (map.getMinInt() < m_min[it_dim])
          m_min[it_dim] = map.getMinInt();
        if (map.getMaxInt() > m_max[it_dim])
          m_max[it_dim] = map.getMaxInt();
      }
      else if (getLayer(layer_index).type == LayerData::DT_IDENT)
      {
        // cout << "recalculateRanges_" << endl;
        const vector<PeptideIdentification> & peptides =
          getLayer(layer_index).peptides;
        for (vector<PeptideIdentification>::const_iterator it =
               peptides.begin(); it != peptides.end(); ++it)
        {
          if (!it->getHits().empty())
          {
            DoubleReal rt = (DoubleReal) it->getMetaValue("RT");
            DoubleReal mz = getIdentificationMZ_(layer_index, *it);
            if (mz < m_min[mz_dim])
              m_min[mz_dim] = mz;
            if (mz > m_max[mz_dim])
              m_max[mz_dim] = mz;
            if (rt < m_min[rt_dim])
              m_min[rt_dim] = rt;
            if (rt > m_max[rt_dim])
              m_max[rt_dim] = rt;
          }
        }
      }
    }
    //Add 1% margin to RT in order to display all the data
    DoubleReal margin = 0.01 * std::max(1.0, m_max[rt_dim] - m_min[rt_dim]);
    m_min[rt_dim] -= margin;
    m_max[rt_dim] += margin;

    overall_data_range_.setMin(m_min);
    overall_data_range_.setMax(m_max);
  }

  DoubleReal SpectrumCanvas::getSnapFactor()
  {
    return snap_factors_[0];
  }

  DoubleReal SpectrumCanvas::getPercentageFactor()
  {
    return percentage_factor_;
  }

  void SpectrumCanvas::recalculateSnapFactor_()
  {

  }

  void SpectrumCanvas::horizontalScrollBarChange(int /*value*/)
  {

  }

  void SpectrumCanvas::verticalScrollBarChange(int /*value*/)
  {

  }

  void SpectrumCanvas::update_(const char *)
  {
    update();
  }

  //this does not work anymore, probably due to Qt::StrongFocus :(
  void SpectrumCanvas::focusOutEvent(QFocusEvent * /*e*/)
  {
    // Alt/Shift pressed and focus lost => change back action mode
    if (action_mode_ != AM_TRANSLATE)
    {
      action_mode_ = AM_TRANSLATE;
      emit actionModeChange();
    }

    //reset peaks
    selected_peak_.clear();
    measurement_start_.clear();

    //update
    update_(__PRETTY_FUNCTION__);
  }

  void SpectrumCanvas::leaveEvent(QEvent * /*e*/)
  {
    //release keyboard, when the mouse pointer leaves
    releaseKeyboard();
  }

  void SpectrumCanvas::enterEvent(QEvent * /*e*/)
  {
    //grab keyboard, as we need to handle key presses
    grabKeyboard();
  }

  void SpectrumCanvas::keyReleaseEvent(QKeyEvent * e)
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

  void SpectrumCanvas::keyPressEvent(QKeyEvent * e)
  {
    // Alt/Shift pressed => change action mode
    if (e->key() == Qt::Key_Control)
    {
      e->accept();
      action_mode_ = AM_ZOOM;
      emit actionModeChange();
    }
    else if (e->key() == Qt::Key_Shift)
    {
      e->accept();
      action_mode_ = AM_MEASURE;
      emit actionModeChange();
    }

    // CTRL+/CTRL- => Zoom stack
    if ((e->modifiers() & Qt::ControlModifier) && (e->key() == Qt::Key_Plus))
    {
      e->accept();
      zoomForward_();
    }
    else if ((e->modifiers() & Qt::ControlModifier) && (e->key() == Qt::Key_Minus))
    {
      e->accept();
      zoomBack_();
    }
    // Arrow keys => translate
    else if (e->key() == Qt::Key_Left)
    {
      e->accept();
      translateLeft_();
    }
    else if (e->key() == Qt::Key_Right)
    {
      e->accept();
      translateRight_();
    }
    else if (e->key() == Qt::Key_Up)
    {
      e->accept();
      translateForward_();
    }
    else if (e->key() == Qt::Key_Down)
    {
      e->accept();
      translateBackward_();
    }
    //Backspace to reset zoom
    else if (e->key() == Qt::Key_Backspace)
    {
      e->accept();
      resetZoom();
    }

    // CTRL+ALT+T => activate timing mode
    if ((e->modifiers() & Qt::ControlModifier) && (e->modifiers() & Qt::AltModifier) && (e->key() == Qt::Key_T))
    {
      e->accept();
      show_timing_ = !show_timing_;
    }

    releaseKeyboard();    // ensure that the key event is passed on to parent widget
    e->ignore();
  }

  void SpectrumCanvas::translateLeft_()
  {
  }

  void SpectrumCanvas::translateRight_()
  {
  }

  void SpectrumCanvas::translateForward_()
  {
  }

  void SpectrumCanvas::translateBackward_()
  {
  }

  void SpectrumCanvas::setAdditionalContextMenu(QMenu * menu)
  {
    context_add_ = menu;
  }

  void SpectrumCanvas::getVisiblePeakData(ExperimentType & map) const
  {
    //clear output experiment
    map.clear(true);

    const LayerData & layer = getCurrentLayer();
    if (layer.type == LayerData::DT_PEAK)
    {
      const AreaType & area = getVisibleArea();
      const ExperimentType & peaks = *layer.getPeakData();
      //copy experimental settings
      map.ExperimentalSettings::operator=(peaks);
      //reserve space for the correct number of spectra in RT range
      ExperimentType::ConstIterator begin = layer.getPeakData()->RTBegin(area.minPosition()[1]);
      ExperimentType::ConstIterator end = layer.getPeakData()->RTEnd(area.maxPosition()[1]);

      //Exception for Spectrum1DCanvas, here we copy the currently visualized spectrum
      bool is_1d = (getName() == "Spectrum1DCanvas");
      if (is_1d)
      {
        begin = layer.getPeakData()->begin() + layer.getCurrentSpectrumIndex();
        end = begin + 1;
      }

      map.reserve(end - begin);
      //copy spectra
      for (ExperimentType::ConstIterator it = begin; it != end; ++it)
      {
        SpectrumType spectrum;
        //copy spectrum meta information
        spectrum.SpectrumSettings::operator=(* it);
        spectrum.setRT(it->getRT());
        spectrum.setMSLevel(it->getMSLevel());
        spectrum.setPrecursors(it->getPrecursors());
        //copy peak information
        if (!is_1d && it->getMSLevel() > 1 && !it->getPrecursors().empty())       //MS^n (n>1) spectra are copied if their precursor is in the m/z range
        {
          if (it->getPrecursors()[0].getMZ() >= area.minPosition()[0] && it->getPrecursors()[0].getMZ() <= area.maxPosition()[0])
          {
            spectrum.insert(spectrum.begin(), it->begin(), it->end());
          }
        }
        else         // MS1(0) spectra are cropped to the m/z range
        {
          for (SpectrumType::ConstIterator it2 = it->MZBegin(area.minPosition()[0]); it2 != it->MZEnd(area.maxPosition()[0]); ++it2)
          {
            if (layer.filters.passes(*it, it2 - it->begin()))
            {
              spectrum.push_back(*it2);
            }
          }
        }
        map.push_back(spectrum);
      }
    }
    else if (layer.type == LayerData::DT_CHROMATOGRAM)
    {
      //TODO CHROM
    }
  }

  void SpectrumCanvas::getVisibleFeatureData(FeatureMapType & map) const
  {
    //clear output experiment
    map.clear(true);

    const LayerData & layer = getCurrentLayer();
    if (layer.type == LayerData::DT_FEATURE)
    {
      //copy meta data
      map.setIdentifier(layer.getFeatureMap()->getIdentifier());
      map.setProteinIdentifications(layer.getFeatureMap()->getProteinIdentifications());
      //Visible area
      DoubleReal min_rt = getVisibleArea().minPosition()[1];
      DoubleReal max_rt = getVisibleArea().maxPosition()[1];
      DoubleReal min_mz = getVisibleArea().minPosition()[0];
      DoubleReal max_mz = getVisibleArea().maxPosition()[0];
      //copy features
      for (FeatureMapType::ConstIterator it = layer.getFeatureMap()->begin(); it != layer.getFeatureMap()->end(); ++it)
      {
        if (layer.filters.passes(*it)
           && it->getRT() >= min_rt
           && it->getRT() <= max_rt
           && it->getMZ() >= min_mz
           && it->getMZ() <= max_mz)
        {
          map.push_back(*it);
        }
      }
    }
  }

  void SpectrumCanvas::getVisibleConsensusData(ConsensusMapType & map) const
  {
    //clear output experiment
    map.clear(true);

    const LayerData & layer = getCurrentLayer();
    if (layer.type == LayerData::DT_CONSENSUS)
    {
      //copy file descriptions
      map.getFileDescriptions() = layer.getConsensusMap()->getFileDescriptions();
      //Visible area
      DoubleReal min_rt = getVisibleArea().minPosition()[1];
      DoubleReal max_rt = getVisibleArea().maxPosition()[1];
      DoubleReal min_mz = getVisibleArea().minPosition()[0];
      DoubleReal max_mz = getVisibleArea().maxPosition()[0];
      //copy features
      for (ConsensusMapType::ConstIterator it = layer.getConsensusMap()->begin(); it != layer.getConsensusMap()->end(); ++it)
      {
        if (layer.filters.passes(*it)
           && it->getRT() >= min_rt
           && it->getRT() <= max_rt
           && it->getMZ() >= min_mz
           && it->getMZ() <= max_mz)
        {
          map.push_back(*it);
        }
      }
    }
  }

  void SpectrumCanvas::getVisibleIdentifications(vector<PeptideIdentification> &
                                                 peptides) const
  {
    //clear output experiment
    peptides.clear();

    const LayerData & layer = getCurrentLayer();
    if (layer.type == LayerData::DT_IDENT)
    {
      //Visible area
      DoubleReal min_rt = getVisibleArea().minPosition()[1];
      DoubleReal max_rt = getVisibleArea().maxPosition()[1];
      DoubleReal min_mz = getVisibleArea().minPosition()[0];
      DoubleReal max_mz = getVisibleArea().maxPosition()[0];
      //copy features
      for (vector<PeptideIdentification>::const_iterator it =
             layer.peptides.begin(); it != layer.peptides.end(); ++it)
      {
        DoubleReal rt = (DoubleReal) it->getMetaValue("RT");
        DoubleReal mz = getIdentificationMZ_(current_layer_, *it);
        // TODO: if (layer.filters.passes(*it) && ...)
        if ((rt >= min_rt) && (rt <= max_rt) &&
            (mz >= min_mz) && (mz <= max_mz))
        {
          peptides.push_back(*it);
        }
      }
    }
  }

  void SpectrumCanvas::showMetaData(bool modifiable, Int index)
  {
    LayerData & layer = getCurrentLayer_();

    MetaDataBrowser dlg(modifiable, this);
    if (index == -1)
    {
      if (layer.type == LayerData::DT_PEAK)
      {
        dlg.add(*layer.getPeakData());
        //Exception for Spectrum1DCanvas, here we add the meta data of the one spectrum
        if (getName() == "Spectrum1DCanvas")
        {
          dlg.add((*layer.getPeakData())[layer.getCurrentSpectrumIndex()]);
        }
      }
      else if (layer.type == LayerData::DT_FEATURE)
      {
        dlg.add(*layer.getFeatureMap());
      }
      else if (layer.type == LayerData::DT_CONSENSUS)
      {
        dlg.add(*layer.getConsensusMap());
      }
      else if (layer.type == LayerData::DT_CHROMATOGRAM)
      {
        //TODO CHROM
      }
      else if (layer.type == LayerData::DT_IDENT)
      {
        // TODO IDENT
      }
    }
    else     //show element meta data
    {
      if (layer.type == LayerData::DT_PEAK)
      {
        dlg.add((*layer.getPeakData())[index]);
      }
      else if (layer.type == LayerData::DT_FEATURE)
      {
        dlg.add((*layer.getFeatureMap())[index]);
      }
      else if (layer.type == LayerData::DT_CONSENSUS)
      {
        dlg.add((*layer.getConsensusMap())[index]);
      }
      else if (layer.type == LayerData::DT_CHROMATOGRAM)
      {
        //TODO CHROM
      }
      else if (layer.type == LayerData::DT_IDENT)
      {
        // TODO IDENT
      }
    }

    //if the meta data was modified, set the flag
    if (modifiable && dlg.exec())
    {
      modificationStatus_(activeLayerIndex(), true);
    }
  }

  void SpectrumCanvas::updateCursor_()
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

  void SpectrumCanvas::modificationStatus_(Size layer_index, bool modified)
  {
    LayerData & layer = getLayer_(layer_index);
    if (layer.modified != modified)
    {
      layer.modified = modified;
      #ifdef DEBUG_TOPPVIEW
      cout << "BEGIN " << __PRETTY_FUNCTION__ << endl;
      cout << "emit: layerModificationChange" << endl;
      cout << "END " << __PRETTY_FUNCTION__ << endl;
      #endif
      emit layerModficationChange(activeLayerIndex(), modified);
    }
  }

  void SpectrumCanvas::drawText_(QPainter & painter, QStringList text)
  {
    painter.save();

    //font
    QFont font("Courier");
    painter.setFont(font);

    //determine width and height of the box we need
    QFontMetrics metrics(painter.font());
    int line_spacing = metrics.lineSpacing();
    int height = 6 + text.size() * line_spacing;
    int width = 4;
    for (int i = 0; i < text.size(); ++i)
    {
      width = std::max(width, 4 + metrics.width(text[i]));
    }

    //draw backgrond for text
    painter.fillRect(2, 3, width, height, QColor(255, 255, 255, 200));

    //draw text
    painter.setPen(Qt::black);
    for (int i = 0; i < text.size(); ++i)
    {
      painter.drawText(3, 3 + (i + 1) * line_spacing, text[i]);
    }

    painter.restore();
  }

  DoubleReal SpectrumCanvas::getIdentificationMZ_(const Size layer_index,
                                                  const PeptideIdentification &
                                                  peptide) const
  {
    if (getLayerFlag(layer_index, LayerData::I_PEPTIDEMZ))
    {
      const PeptideHit & hit = peptide.getHits().front();
      Int charge = hit.getCharge();
      return hit.getSequence().getMonoWeight(Residue::Full, charge) / charge;
    }
    else
    {
      return (DoubleReal) peptide.getMetaValue("MZ");
    }
  }

} //namespace
