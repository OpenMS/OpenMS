// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
    spectrum_widget_(nullptr),
    percentage_factor_(1.0),
    snap_factors_(1, 1.0),
    rubber_band_(QRubberBand::Rectangle, this),
    context_add_(nullptr),
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

  void SpectrumCanvas::setFilters(const DataFilters & filters)
  {
    //set filters
    layers_[current_layer_].filters = filters;
    //update the content
    update_buffer_ = true;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void SpectrumCanvas::showGridLines(bool show)
  {
    show_grid_ = show;
    update_buffer_ = true;
    update_(OPENMS_PRETTY_FUNCTION);
  }

  void SpectrumCanvas::intensityModeChange_()
  {
    recalculateSnapFactor_();
    update_buffer_ = true;
    update_(OPENMS_PRETTY_FUNCTION);
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
    update_(OPENMS_PRETTY_FUNCTION);
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
      update_(OPENMS_PRETTY_FUNCTION);
    }
  }

  void SpectrumCanvas::updateScrollbars_()
  {

  }

  void SpectrumCanvas::wheelEvent(QWheelEvent * e)
  {
    zoom_(e->x(), e->y(), e->delta() > 0);
    e->accept();
  }

  void SpectrumCanvas::zoom_(int x, int y, bool zoom_in)
  {
    const PointType::CoordinateType zoom_factor = zoom_in ? 0.8 : 1.0 / 0.8;
    AreaType new_area;
    for (int dim = 0; dim < AreaType::DIMENSION; dim++)
    {
      // don't assign to "new_area.min_"/"max_" immediately, as this can lead to strange crashes at the min/max calls below on some platforms
      // (GCC 4.6.3; faulty out-of-order execution?):
      AreaType::CoordinateType coef = ((dim == 0) == isMzToXAxis()) ? (AreaType::CoordinateType(x) / width()) : (AreaType::CoordinateType(height() - y) / height());
      AreaType::CoordinateType min_pos = visible_area_.min_[dim] + (1.0 - zoom_factor) * (visible_area_.max_[dim] - visible_area_.min_[dim]) * coef;
      AreaType::CoordinateType max_pos = min_pos + zoom_factor * (visible_area_.max_[dim] - visible_area_.min_[dim]);
      new_area.min_[dim] = max(min_pos, overall_data_range_.min_[dim]);
      new_area.max_[dim] = min(max_pos, overall_data_range_.max_[dim]);
    }
    if (new_area != visible_area_)
    {
      zoomAdd_(new_area);
      zoom_pos_ = --zoom_stack_.end(); // set to last position
      changeVisibleArea_(*zoom_pos_);
    }
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
    //cout << OPENMS_PRETTY_FUNCTION << endl;
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

    if (layers_.back().getPeakData()->getChromatograms().size() != 0 
        && layers_.back().getPeakData()->size() != 0)
    {
      // TODO : handle this case better
      LOG_WARN << "Your input data contains chromatograms and spectra, falling back to display spectra only." << std::endl;
    }

    if (layers_.back().getPeakData()->getChromatograms().size() != 0 
        && layers_.back().getPeakData()->size() == 0)
    {
      layers_.back().type = LayerData::DT_CHROMATOGRAM;
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
    OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::setLayerName(i, name) index overflow");
    getLayer_(i).name = name;
    if (i == 0 && spectrum_widget_) 
    {
      spectrum_widget_->setWindowTitle(name.toQString());
    }
  }

  String SpectrumCanvas::getLayerName(const Size i)
  {
    OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::getLayerName(i) index overflow");
    return getLayer_(i).name;
  }

  void SpectrumCanvas::changeVisibility(Size i, bool b)
  {
    OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::changeVisibility(i, b) index overflow");
    LayerData & layer = getLayer_(i);
    if (layer.visible != b)
    {
      layer.visible = b;
      update_buffer_ = true;
      update_(OPENMS_PRETTY_FUNCTION);
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
      update_(OPENMS_PRETTY_FUNCTION);
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
        if (map.getMinMZ() < m_min[mz_dim]) m_min[mz_dim] = map.getMinMZ();
        if (map.getMaxMZ() > m_max[mz_dim]) m_max[mz_dim] = map.getMaxMZ();
        if (map.getMinRT() < m_min[rt_dim]) m_min[rt_dim] = map.getMinRT();
        if (map.getMaxRT() > m_max[rt_dim]) m_max[rt_dim] = map.getMaxRT();
        if (map.getMinInt() < m_min[it_dim]) m_min[it_dim] = map.getMinInt();
        if (map.getMaxInt() > m_max[it_dim]) m_max[it_dim] = map.getMaxInt();
      }
      else if (getLayer(layer_index).type == LayerData::DT_FEATURE)
      {
        const FeatureMapType & map = *getLayer(layer_index).getFeatureMap();
        if (map.getMin()[1] < m_min[mz_dim]) m_min[mz_dim] = map.getMin()[1];
        if (map.getMax()[1] > m_max[mz_dim]) m_max[mz_dim] = map.getMax()[1];
        if (map.getMin()[0] < m_min[rt_dim]) m_min[rt_dim] = map.getMin()[0];
        if (map.getMax()[0] > m_max[rt_dim]) m_max[rt_dim] = map.getMax()[0];
        if (map.getMinInt() < m_min[it_dim]) m_min[it_dim] = map.getMinInt();
        if (map.getMaxInt() > m_max[it_dim]) m_max[it_dim] = map.getMaxInt();
      }
      else if (getLayer(layer_index).type == LayerData::DT_CONSENSUS)
      {
        const ConsensusMapType & map = *getLayer(layer_index).getConsensusMap();
        if (map.getMin()[1] < m_min[mz_dim]) m_min[mz_dim] = map.getMin()[1];
        if (map.getMax()[1] > m_max[mz_dim]) m_max[mz_dim] = map.getMax()[1];
        if (map.getMin()[0] < m_min[rt_dim]) m_min[rt_dim] = map.getMin()[0];
        if (map.getMax()[0] > m_max[rt_dim]) m_max[rt_dim] = map.getMax()[0];
        if (map.getMinInt() < m_min[it_dim]) m_min[it_dim] = map.getMinInt();
        if (map.getMaxInt() > m_max[it_dim]) m_max[it_dim] = map.getMaxInt();
      }
      else if (getLayer(layer_index).type == LayerData::DT_IDENT)
      {
        const vector<PeptideIdentification> & peptides =
          getLayer(layer_index).peptides;
        for (vector<PeptideIdentification>::const_iterator it =
               peptides.begin(); it != peptides.end(); ++it)
        {
          double rt = it->getRT();
          double mz = getIdentificationMZ_(layer_index, *it);
          if (mz < m_min[mz_dim]) m_min[mz_dim] = mz;
          if (mz > m_max[mz_dim]) m_max[mz_dim] = mz;
          if (rt < m_min[rt_dim]) m_min[rt_dim] = rt;
          if (rt > m_max[rt_dim]) m_max[rt_dim] = rt;
        }
      }
    }
    overall_data_range_.setMin(m_min);
    overall_data_range_.setMax(m_max);

    // add 4% margin (2% left, 2% right) to RT, m/z and intensity
    overall_data_range_.extend(1.04);

    // set minimum intensity to 0
    DRange<3>::PositionType new_min = overall_data_range_.minPosition();
    new_min[it_dim] = 0;
    overall_data_range_.setMin(new_min);
  }

  double SpectrumCanvas::getSnapFactor()
  {
    return snap_factors_[0];
  }

  double SpectrumCanvas::getPercentageFactor()
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
    update_(OPENMS_PRETTY_FUNCTION);
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
    else if ((e->modifiers() & Qt::ControlModifier) && (e->key() == Qt::Key_Plus)) // do not use (e->modifiers() == Qt::ControlModifier) to target Ctrl exclusively, since +/- might(!) also trigger the Qt::KeypadModifier
    { // CTRL+Plus => Zoom stack
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

  void SpectrumCanvas::translateLeft_(Qt::KeyboardModifiers /*m*/)
  {
  }

  void SpectrumCanvas::translateRight_(Qt::KeyboardModifiers /*m*/)
  {
  }

  void SpectrumCanvas::translateForward_()
  {
  }

  void SpectrumCanvas::translateBackward_()
  {
  }

  void SpectrumCanvas::setAdditionalContextMenu(QMenu* menu)
  {
    context_add_ = menu;
  }

  void SpectrumCanvas::getVisiblePeakData(ExperimentType& map) const
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
            map.addSpectrum(spectrum);
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
          map.addSpectrum(spectrum);
        }
        // do not use map.addSpectrum() here, otherwise empty spectra which did not pass the filters above will be added
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
      double min_rt = getVisibleArea().minPosition()[1];
      double max_rt = getVisibleArea().maxPosition()[1];
      double min_mz = getVisibleArea().minPosition()[0];
      double max_mz = getVisibleArea().maxPosition()[0];
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
      double min_rt = getVisibleArea().minPosition()[1];
      double max_rt = getVisibleArea().maxPosition()[1];
      double min_mz = getVisibleArea().minPosition()[0];
      double max_mz = getVisibleArea().maxPosition()[0];
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
      double min_rt = getVisibleArea().minPosition()[1];
      double max_rt = getVisibleArea().maxPosition()[1];
      double min_mz = getVisibleArea().minPosition()[0];
      double max_mz = getVisibleArea().maxPosition()[0];
      //copy features
      for (vector<PeptideIdentification>::const_iterator it =
             layer.peptides.begin(); it != layer.peptides.end(); ++it)
      {
        double rt = it->getRT();
        double mz = getIdentificationMZ_(current_layer_, *it);
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
      cout << "BEGIN " << OPENMS_PRETTY_FUNCTION << endl;
      cout << "emit: layerModificationChange" << endl;
      cout << "END " << OPENMS_PRETTY_FUNCTION << endl;
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
      width = max(width, 4 + metrics.width(text[i]));
    }

    //draw background for text
    painter.fillRect(2, 3, width, height, QColor(255, 255, 255, 200));

    //draw text
    painter.setPen(Qt::black);
    for (int i = 0; i < text.size(); ++i)
    {
      painter.drawText(3, 3 + (i + 1) * line_spacing, text[i]);
    }
    painter.restore();
  }

  double SpectrumCanvas::getIdentificationMZ_(const Size layer_index,
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
      return peptide.getMZ();
    }
  }

} //namespace
