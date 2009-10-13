// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>

#include <boost/math/special_functions/fpclassify.hpp>

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
	SpectrumCanvas::SpectrumCanvas(const Param& /*preferences*/, QWidget* parent)
		: QWidget(parent),
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
			snap_factors_(1,1.0),
			rubber_band_(QRubberBand::Rectangle,this),
			watcher_(0),
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
			
		setMinimumSize(200,200);
		setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
	  
	  //reserve enough space to avoid copying layer data
	  layers_.reserve(10);
	  
	  //set common defaults for all canvases
    defaults_.setValue("default_path", ".", "Default path for loading/storing data.");
    defaults_.setValue("on_file_change", "ask", "What action to take, when a data file changes. Do nothing, update automatically or ask the user.");
    defaults_.setValidStrings("on_file_change",StringList::create("none,ask,update automatically"));
    
    //create file system watcher
    watcher_ = new FileWatcher(this);
    connect(watcher_,SIGNAL(fileChanged(const String&)),this,SLOT(fileChanged_(const String&)));
    
    //Set focus policy in order to get keyboard events
	  
	  
	  //Set 'whats this' text
	  setWhatsThis("Translate: Translate mode is activated by default. Hold down the left mouse key and move the mouse to translate. Arrow keys can be used for translation independent of the current mode.\n\n"
	  						 "Zoom: Zoom mode is activated with the CTRL key. CTRL+/CTRL- are used to traverse the zoom stack (or mouse wheel). Pressing Backspace resets the zoom.\n\n"
	  						 "Measure: Measure mode is activated with the SHIFT key. To measure the distace between data points, press the left mouse button on a point and drag the mouse to another point.\n\n"
								 );
		
		//set move cursor and connect signal that updates the cursor automatically
		updateCursor_();
		connect(this,SIGNAL(actionModeChange()),this,SLOT(updateCursor_()));
	}

	SpectrumCanvas::~SpectrumCanvas()
	{
		//cout << "DEST SpectrumCanvas" << endl;
	}

	void SpectrumCanvas::resizeEvent(QResizeEvent* /* e */)
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

	void SpectrumCanvas::setFilters(const DataFilters& filters)
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
	
	void SpectrumCanvas::changeVisibleArea_(const AreaType& new_area, bool repaint, bool add_to_stack)
	{
		//store old zoom state
		if (add_to_stack)
		{
			// if we scrolled in between zooming we want to store the last position before zooming as well
			if (	 (zoom_stack_.size()>0)
					&& (zoom_stack_.back()!=visible_area_))
			{
				zoomAdd_(visible_area_);
			}
			// add current zoom
			zoomAdd_(new_area);
		}
		
		if (new_area!=visible_area_)
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

	void SpectrumCanvas::wheelEvent(QWheelEvent* e)
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
		if (zoom_pos_!=zoom_stack_.begin())
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
		if (zoom_pos_!=zoom_stack_.end() && (zoom_pos_+1)!=zoom_stack_.end())
		{
			++zoom_pos_;
			changeVisibleArea_(*zoom_pos_);
		}
		//cout << " - pos after:" << (zoom_pos_-zoom_stack_.begin()) << endl;
	}
	
	void SpectrumCanvas::zoomAdd_(const AreaType& area)
	{
		//cout << "Adding to stack" << endl;
		//cout << " - pos before:" << (zoom_pos_-zoom_stack_.begin()) << endl;
		//cout << " - size before:" << zoom_stack_.size() <<endl;
		if (zoom_pos_!=zoom_stack_.end() && (zoom_pos_+1)!=zoom_stack_.end())
		{
			//cout << " - removing from:" << ((zoom_pos_+1)-zoom_stack_.begin()) << endl;
			zoom_stack_.erase(zoom_pos_+1,zoom_stack_.end());
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
		changeVisibleArea_(tmp,repaint,true);
	}
	
	void SpectrumCanvas::setVisibleArea(AreaType area)
	{
		//cout << __PRETTY_FUNCTION__ << endl;
		changeVisibleArea_(area);
	}
	
	
	void SpectrumCanvas::paintGridLines_(QPainter& painter)
	{	
		if (!show_grid_ || !spectrum_widget_) return;

		QPen p1(QColor(130,130,130));
		p1.setStyle(Qt::DashLine);
		QPen p2(QColor(170,170,170));
		p2.setStyle(Qt::DashLine);
		QPen p3(QColor(230,230,230));
		p3.setStyle(Qt::DashLine);
	
		painter.save();

		unsigned int xl, xh, yl, yh; //width/height of the diagram area, x, y coordinates of lo/hi x,y values
	
		xl = 0;
		xh = width();

		yl = height();
		yh = 0;
	
		// drawing of grid lines and associated text	
		for (Size j = 0; j != spectrum_widget_->xAxis()->gridLines().size() ; j++) 
		{
			// style definitions
			switch(j)
			{
				case 0:	// style settings for big intervals 
					painter.setPen(p1);
					break;
				case 1:	// style settings for small intervals
					painter.setPen(p2);
					break;
				case 2: // style settings for smalles intervals
					painter.setPen(p3);
					break;
				default:
					std::cout << "empty vertical grid line vector error!" << std::endl;
					painter.setPen(QPen(QColor(0,0,0)));
					break;
			}

			int x;
			for (std::vector<double>::const_iterator it = spectrum_widget_->xAxis()->gridLines()[j].begin(); it != spectrum_widget_->xAxis()->gridLines()[j].end(); it++) 
			{
				x = static_cast<int>(Math::intervalTransformation(*it, spectrum_widget_->xAxis()->getAxisMinimum(), spectrum_widget_->xAxis()->getAxisMaximum(), xl, xh));
				painter.drawLine(x, yl, x, yh);
			}
		}
		
		for (Size j = 0; j != spectrum_widget_->yAxis()->gridLines().size() ; j++) 
		{

			// style definitions
			switch(j)
			{
				case 0:	// style settings for big intervals 
					painter.setPen(p1);
					break;
				case 1:	// style settings for small intervals
					painter.setPen(p2);
					break;
				case 2: // style settings for smalles intervals
					painter.setPen(p3);
					break;
				default:
					std::cout << "empty vertical grid line vector error!" << std::endl;
					painter.setPen(QPen(QColor(0,0,0)));
					break;
			}

			int y;
			for (std::vector<double>::const_iterator it = spectrum_widget_->yAxis()->gridLines()[j].begin(); it != spectrum_widget_->yAxis()->gridLines()[j].end(); it++) 
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

	bool SpectrumCanvas::addLayer(ExperimentType& map, const String& filename)
	{	
		layers_.resize(layers_.size()+1);
		layers_.back().param = param_;
		layers_.back().filename = filename;
		layers_.back().peaks.swap(map);
		if (layers_.back().peaks.getChromatograms().size()!=0)
		{
			Size num_chrom(0);
			for (Size i = 0; i != layers_.back().peaks.getChromatograms().size(); ++i)
			{
				if (layers_.back().peaks.getChromatograms()[i].getChromatogramType() == ChromatogramSettings::SELECTED_ION_CURRENT_CHROMATOGRAM ||
						layers_.back().peaks.getChromatograms()[i].getChromatogramType() == ChromatogramSettings::SELECTED_REACTION_MONITORING_CHROMATOGRAM)
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

	bool SpectrumCanvas::addLayer(FeatureMapType& map, const String& filename)
	{
		layers_.resize(layers_.size()+1);
		layers_.back().param = param_;
		layers_.back().filename = filename;
		layers_.back().features.swap(map);
		layers_.back().type = LayerData::DT_FEATURE;

		return finishAdding_();
	}

	bool SpectrumCanvas::addLayer(ConsensusMapType& map, const String& filename)
	{
		layers_.resize(layers_.size()+1);
		layers_.back().param = param_;
		layers_.back().filename = filename;
		layers_.back().consensus.swap(map);
		layers_.back().type = LayerData::DT_CONSENSUS;

		return finishAdding_();
	}

	void SpectrumCanvas::setLayerName(Size i, const String& name)
	{ 
		OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::setLayerName(i,name) index overflow");
	  getLayer_(i).name = name; 
		if (i==0 && spectrum_widget_) spectrum_widget_->setWindowTitle(name.toQString());
	}

	void SpectrumCanvas::changeVisibility(Size i, bool b)
	{
		OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::changeVisibility(i,b) index overflow");
		LayerData& layer = getLayer_(i);
		if (layer.visible!=b)
		{
			layer.visible=b;
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
	}

	void SpectrumCanvas::changeLayerFilterState(Size i, bool b)
	{
		OPENMS_PRECONDITION(i < layers_.size(), "SpectrumCanvas::changeVisibility(i,b) index overflow");
		LayerData& layer = getLayer_(i);
		if (layer.filters.isActive()!=b)
		{
			layer.filters.setActive(b);
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
	}

  const DRange<3>& SpectrumCanvas::getDataRange()
  {
  	return overall_data_range_;
  }
	
	void SpectrumCanvas::recalculateRanges_(UInt mz_dim, UInt rt_dim, UInt it_dim)
	{
		overall_data_range_ = DRange<3>::empty;
		DRange<3>::PositionType min = overall_data_range_.min();
		DRange<3>::PositionType max = overall_data_range_.max();
		
		for (Size layer_index=0; layer_index< getLayerCount(); ++layer_index)
		{
			if (getLayer(layer_index).type==LayerData::DT_PEAK || getLayer(layer_index).type==LayerData::DT_CHROMATOGRAM)
			{
				const ExperimentType& map = getLayer(layer_index).peaks;
				if (map.getMinMZ() < min[mz_dim]) min[mz_dim] = map.getMinMZ();
				if (map.getMaxMZ() > max[mz_dim]) max[mz_dim] = map.getMaxMZ();
				if (map.getMinRT() < min[rt_dim]) min[rt_dim] = map.getMinRT();
				if (map.getMaxRT() > max[rt_dim]) max[rt_dim] = map.getMaxRT();
				if (map.getMinInt() < min[it_dim]) min[it_dim] = map.getMinInt();
				if (map.getMaxInt() > max[it_dim]) max[it_dim] = map.getMaxInt();
			}
			else if (getLayer(layer_index).type==LayerData::DT_FEATURE)
			{
				const FeatureMapType& map = getLayer(layer_index).features;
				if (map.getMin()[1] < min[mz_dim]) min[mz_dim] = map.getMin()[1];
				if (map.getMax()[1] > max[mz_dim]) max[mz_dim] = map.getMax()[1];
				if (map.getMin()[0] < min[rt_dim]) min[rt_dim] = map.getMin()[0];
				if (map.getMax()[0] > max[rt_dim]) max[rt_dim] = map.getMax()[0];
				if (map.getMinInt() < min[it_dim]) min[it_dim] = map.getMinInt();
				if (map.getMaxInt() > max[it_dim]) max[it_dim] = map.getMaxInt();
			}
			else
			{
				const ConsensusMapType& map = getLayer(layer_index).consensus;
				if (map.getMin()[1] < min[mz_dim]) min[mz_dim] = map.getMin()[1];
				if (map.getMax()[1] > max[mz_dim]) max[mz_dim] = map.getMax()[1];
				if (map.getMin()[0] < min[rt_dim]) min[rt_dim] = map.getMin()[0];
				if (map.getMax()[0] > max[rt_dim]) max[rt_dim] = map.getMax()[0];
				if (map.getMinInt() < min[it_dim]) min[it_dim] = map.getMinInt();
				if (map.getMaxInt() > max[it_dim]) max[it_dim] = map.getMaxInt();
			}
		}
		//Add 1% margin to RT in order to display all the data
		DoubleReal margin = 0.01*std::max(1.0, max[rt_dim] - min[rt_dim]);
		min[rt_dim] -= margin;
		max[rt_dim] += margin;
		//Add 1% margin to MZ in order to display all the data
		margin = 0.01*std::max(1.0, max[mz_dim] - min[mz_dim]);
		min[mz_dim] -= margin;
		max[mz_dim] += margin;
		
		overall_data_range_.setMin(min);
		overall_data_range_.setMax(max);
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
	
	void SpectrumCanvas::update_(const char*
#ifdef DEBUG_UPDATE_
			caller_name)
	{
		cout << "Spectrum3DCanvas::update_ from '" << caller_name << "'" << endl;
#else
		)
	{
#endif
		update();
	}
	
	void SpectrumCanvas::fileChanged_(const String& filename)
  {
		//look up all layers that contain data of the file
		UInt updatable_layers = 0;
		for (UInt j=0; j<getLayerCount(); ++j)
		{	
			//cout << "  Layer: " << j << " " << getLayer(j).filename << endl;
			if (getLayer(j).filename == filename)
			{
				++updatable_layers;
				bool update = false;
				if ((String)(param_.getValue("on_file_change"))=="update automatically") //automatically update
				{
					update = true;
				}
				else if ((String)(param_.getValue("on_file_change"))=="ask") //ask the user if the layer should be updated
				{
					QMessageBox msg_box;
					QAbstractButton* ok = msg_box.addButton(QMessageBox::Ok);
					msg_box.addButton(QMessageBox::Cancel);
					msg_box.setWindowTitle("Layer data changed");
					msg_box.setText((String("The data file of layer '") + getLayer(j).filename + "' has changed.<BR>Update the layer?").toQString());
					msg_box.exec();
					if (msg_box.clickedButton() == ok)
					{
						update = true;
					}
				}
				//update the layer if the user choosed to do so
				if (update)
				{
					emit sendStatusMessage(String("Updating layer '") + getLayer(j).name + "' (file changed).",0);
					updateLayer_(j);
					emit sendStatusMessage(String("Finished updating layer '") + getLayer(j).name + "'.",5000);
				}
			}
		}
		//remove watchers that are not needed anymore
  	if (updatable_layers==0)
  	{
  		watcher_->removeFile(filename);
  	}  			
	}
	
	//this does not work anymore, probably due to Qt::StrongFocus :(
	void SpectrumCanvas::focusOutEvent(QFocusEvent* /*e*/)
	{
		// Alt/Shift pressed and focus lost => change back action mode
		if (action_mode_!=AM_TRANSLATE)
		{
			action_mode_ = AM_TRANSLATE;
			emit actionModeChange();
		}
		
		//reset peaks
		selected_peak_.clear();
		measurement_start_.clear();
		
		//update
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}

	void SpectrumCanvas::leaveEvent(QEvent* /*e*/)
	{
		//release keyboard, when the mouse pointer leaves
		releaseKeyboard();
	}

	void SpectrumCanvas::enterEvent(QEvent* /*e*/)
	{
		//grab keyboard, as we need to handle key presses
		grabKeyboard();
	}

	void SpectrumCanvas::keyReleaseEvent(QKeyEvent* e)
	{
		// Alt/Shift released => change back action mode
		if (e->key()==Qt::Key_Control || e->key()==Qt::Key_Shift)
		{
			action_mode_ = AM_TRANSLATE;
			emit actionModeChange();
			e->accept();
		}

		e->ignore();
	}

	void SpectrumCanvas::keyPressEvent(QKeyEvent* e)
	{
		// Alt/Shift pressed => change action mode
		if (e->key()==Qt::Key_Control)
		{
			e->accept();
			action_mode_ = AM_ZOOM;
			emit actionModeChange();
		}
		else if (e->key()==Qt::Key_Shift)
		{
			e->accept();
			action_mode_ = AM_MEASURE;
			emit actionModeChange();
		}
		
		// CTRL+/CTRL- => Zoom stack
		if ((e->modifiers() & Qt::ControlModifier) && (e->key()==Qt::Key_Plus))
		{
			e->accept();
			zoomForward_();
		}
		else if ((e->modifiers() & Qt::ControlModifier) && (e->key()==Qt::Key_Minus))
		{
			e->accept();
			zoomBack_();
		}
		
		// Arrow keys => translate
		else if (e->key()==Qt::Key_Left)
		{
			e->accept();
			translateLeft_();
		}
		else if (e->key()==Qt::Key_Right)
		{
			e->accept();
			translateRight_();
		}
		else if (e->key()==Qt::Key_Up)
		{
			e->accept();
			translateForward_();
		}
		else if (e->key()==Qt::Key_Down)
		{
			e->accept();
			translateBackward_();
		}
		
		//Backspace to reset zoom
		else if (e->key()==Qt::Key_Backspace)
		{
			e->accept();
			resetZoom();
		}

		// CTRL+ALT+T => activate timing mode
		if ((e->modifiers() & Qt::ControlModifier) && (e->modifiers() & Qt::AltModifier) && (e->key()==Qt::Key_T))
		{
			e->accept();
			show_timing_ = !show_timing_;
		}
		
		releaseKeyboard();// ensure that the key event is passed on to parent widget
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

	void SpectrumCanvas::setAdditionalContextMenu(QMenu* menu)
	{
	  context_add_ = menu;
	}
	
	void SpectrumCanvas::getVisiblePeakData(ExperimentType& map) const
	{		
		//clear output experiment
		map.clear(true);
		
    const LayerData& layer = getCurrentLayer();
  	if (layer.type==LayerData::DT_PEAK)
  	{
			const AreaType& area = getVisibleArea();
			const ExperimentType& peaks = layer.peaks;
			//copy experimental settings
			map.ExperimentalSettings::operator=(peaks);
			//reserve space for the correct number of spectra in RT range
			ExperimentType::ConstIterator begin = layer.peaks.RTBegin(area.min()[1]);
			ExperimentType::ConstIterator end = layer.peaks.RTEnd(area.max()[1]);
			
			//Exception for Spectrum1DCanvas, here we copy the currently visualized spectrum
			bool is_1d = (getName()=="Spectrum1DCanvas");
			if (is_1d)
			{
				begin = layer.peaks.begin() + layer.current_spectrum;
				end = begin+1;
			}

			map.reserve(end-begin);
			//copy spectra
  		for (ExperimentType::ConstIterator it=begin; it!=end; ++it)
  		{
  			SpectrumType spectrum;
				//copy spectrum meta information
				spectrum.SpectrumSettings::operator=(*it);
				spectrum.setRT(it->getRT());
				spectrum.setMSLevel(it->getMSLevel());
				spectrum.setPrecursors(it->getPrecursors());
				//copy peak information
				if (!is_1d && it->getMSLevel()>1 && !it->getPrecursors().empty()) //MS^n (n>1) spectra are copied if their precursor is in the m/z range
				{
					if (it->getPrecursors()[0].getMZ()>=area.min()[0] && it->getPrecursors()[0].getMZ()<= area.max()[0])
					{
						spectrum.insert(spectrum.begin(), it->begin(), it->end());
					}
				}
				else // MS1(0) spectra are cropped to the m/z range
				{
					for (SpectrumType::ConstIterator it2 = it->MZBegin(area.min()[0]); it2!= it->MZEnd(area.max()[0]); ++it2)
					{
						if (layer.filters.passes(*it,it2-it->begin()))
						{
							spectrum.push_back(*it2);
						}
					}
				}
				map.push_back(spectrum);
  		}
		}
		else if (layer.type==LayerData::DT_CHROMATOGRAM)
		{
			//TODO CHROM
		}
	}

	void SpectrumCanvas::getVisibleFeatureData(FeatureMapType& map) const
	{		
		//clear output experiment
		map.clear(true);
		
    const LayerData& layer = getCurrentLayer();
  	if (layer.type==LayerData::DT_FEATURE)
  	{
			//copy meta data
			map.setIdentifier(layer.features.getIdentifier());
			map.setProteinIdentifications(layer.features.getProteinIdentifications());
			//Visible area
			DoubleReal min_rt = getVisibleArea().min()[1];
			DoubleReal max_rt = getVisibleArea().max()[1];
			DoubleReal min_mz = getVisibleArea().min()[0];
			DoubleReal max_mz = getVisibleArea().max()[0];
			//copy features
  		for (FeatureMapType::ConstIterator it=layer.features.begin(); it!=layer.features.end(); ++it)
  		{
				if ( layer.filters.passes(*it) 
					&& it->getRT() >= min_rt 
					&& it->getRT() <= max_rt 
					&& it->getMZ() >= min_mz 
					&& it->getMZ() <= max_mz )
				{
					map.push_back(*it);
				}
			}
		}
	}

	void SpectrumCanvas::getVisibleConsensusData(ConsensusMapType& map) const
	{		
		//clear output experiment
		map.clear(true);
		
    const LayerData& layer = getCurrentLayer();
  	if (layer.type==LayerData::DT_CONSENSUS)
  	{
			//copy file descriptions
			map.getFileDescriptions() = layer.consensus.getFileDescriptions();
			//Visible area
			DoubleReal min_rt = getVisibleArea().min()[1];
			DoubleReal max_rt = getVisibleArea().max()[1];
			DoubleReal min_mz = getVisibleArea().min()[0];
			DoubleReal max_mz = getVisibleArea().max()[0];
			//copy features
  		for (ConsensusMapType::ConstIterator it=layer.consensus.begin(); it!=layer.consensus.end(); ++it)
  		{
				if ( layer.filters.passes(*it)
					&& it->getRT() >= min_rt 
					&& it->getRT() <= max_rt 
					&& it->getMZ() >= min_mz 
					&& it->getMZ() <= max_mz )
				{
					map.push_back(*it);
				}
			}
		}
	}

	void SpectrumCanvas::showMetaData(bool modifiable, Int index)
  {
		LayerData& layer = getCurrentLayer_();
		
		MetaDataBrowser dlg(modifiable, this);
		if (index==-1)
		{
			if (layer.type==LayerData::DT_PEAK)
			{
				dlg.add(layer.peaks);
				//Exception for Spectrum1DCanvas, here we add the meta data of the one spectrum
				if (getName()=="Spectrum1DCanvas")
				{
					dlg.add(layer.peaks[layer.current_spectrum]);
				}
			}
			else if (layer.type==LayerData::DT_FEATURE)
			{
				dlg.add(layer.features);
			}
			else if (layer.type==LayerData::DT_CONSENSUS)
			{
				dlg.add(layer.consensus);
			}
			else if (layer.type==LayerData::DT_CHROMATOGRAM)
			{
				//TODO CHROM
			}
		}
		else //show element meta data
		{
			if (layer.type==LayerData::DT_PEAK)
			{
					dlg.add(layer.peaks[index]);
			}
			else if (layer.type==LayerData::DT_FEATURE)
			{
				dlg.add(layer.features[index]);
			}
			else if (layer.type==LayerData::DT_CONSENSUS)
			{
				dlg.add(layer.consensus[index]);
			}
			else if (layer.type==LayerData::DT_CHROMATOGRAM)
			{
				//TODO CHROM
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
		switch(action_mode_)
		{
			case AM_TRANSLATE:
				setCursor(QCursor(QPixmap(":/cursor_move.png"),0,0));
				break;
			case AM_ZOOM:
				setCursor(QCursor(QPixmap(":/cursor_zoom.png"),0,0));
				break;
			case AM_MEASURE:
				setCursor(QCursor(QPixmap(":/cursor_measure.png"),0,0));
				break;
		}
	}

	void SpectrumCanvas::modificationStatus_(Size layer_index, bool modified)
	{
		LayerData& layer = getLayer_(layer_index);
		if (layer.modified!=modified)
		{
			layer.modified = modified;
			emit layerModficationChange(activeLayerIndex(), modified);
		}
	}


	void SpectrumCanvas::drawCoordinates_(QPainter& painter, const PeakIndex& peak, bool print_rt)
	{
		if (!peak.isValid()) return;
		
		//determine coordinates;
		DoubleReal mz = 0.0;
		DoubleReal rt = 0.0;
		Real it = 0.0;
		Int charge = 0;
		DoubleReal quality = 0.0;
		if (getCurrentLayer().type==LayerData::DT_FEATURE)
		{
			mz = peak.getFeature(getCurrentLayer().features).getMZ();
			rt = peak.getFeature(getCurrentLayer().features).getRT();
			it = peak.getFeature(getCurrentLayer().features).getIntensity();
			charge  = peak.getFeature(getCurrentLayer().features).getCharge();
			quality = peak.getFeature(getCurrentLayer().features).getOverallQuality();
		}
		else if (getCurrentLayer().type==LayerData::DT_PEAK)
		{
			mz = peak.getPeak(getCurrentLayer().peaks).getMZ();
			rt = peak.getSpectrum(getCurrentLayer().peaks).getRT();
			it = peak.getPeak(getCurrentLayer().peaks).getIntensity();
		}
		else if (getCurrentLayer().type==LayerData::DT_CONSENSUS)
		{
			mz = peak.getFeature(getCurrentLayer().consensus).getMZ();
			rt = peak.getFeature(getCurrentLayer().consensus).getRT();
			it = peak.getFeature(getCurrentLayer().consensus).getIntensity();
			charge  = peak.getFeature(getCurrentLayer().consensus).getCharge();
			quality = peak.getFeature(getCurrentLayer().consensus).getQuality();
		}
		else if (getCurrentLayer().type==LayerData::DT_CHROMATOGRAM)
		{
			//TODO CHROM
		}
		
		//draw text			
		QStringList lines;
		if (print_rt) lines.push_back("RT : " + QString::number(rt,'f',2));
		lines.push_back("m/z: " + QString::number(mz,'f',6));
		lines.push_back("Int: " + QString::number(it,'f',2));
		if (getCurrentLayer().type==LayerData::DT_FEATURE || getCurrentLayer().type==LayerData::DT_CONSENSUS)
		{
			lines.push_back("Charge: " + QString::number(charge));
			lines.push_back("Quality: " + QString::number(quality,'f',4));
		}
		drawText_(painter, lines);
	}

	void SpectrumCanvas::drawDeltas_(QPainter& painter, const PeakIndex& start, const PeakIndex& end, bool print_rt)
	{
		if (!start.isValid()) return;
		
		//determine coordinates;
		DoubleReal mz = 0.0;
		DoubleReal rt = 0.0;
		Real it = 0.0;
		if (getCurrentLayer().type==LayerData::DT_FEATURE)
		{
			if (end.isValid())
			{
				mz = end.getFeature(getCurrentLayer().features).getMZ() - start.getFeature(getCurrentLayer().features).getMZ();
				rt = end.getFeature(getCurrentLayer().features).getRT() - start.getFeature(getCurrentLayer().features).getRT();
				it = end.getFeature(getCurrentLayer().features).getIntensity() / start.getFeature(getCurrentLayer().features).getIntensity();
			}
			else
			{
				PointType point = widgetToData_(last_mouse_pos_);
				mz = point[0] - start.getFeature(getCurrentLayer().features).getMZ();
				rt = point[1] - start.getFeature(getCurrentLayer().features).getRT();
				it = std::numeric_limits<DoubleReal>::quiet_NaN();
			}
		}
		else if (getCurrentLayer().type==LayerData::DT_PEAK)
		{
			if (end.isValid())
			{
				mz = end.getPeak(getCurrentLayer().peaks).getMZ() - start.getPeak(getCurrentLayer().peaks).getMZ();
				rt = end.getSpectrum(getCurrentLayer().peaks).getRT() - start.getSpectrum(getCurrentLayer().peaks).getRT();
				it = end.getPeak(getCurrentLayer().peaks).getIntensity() / start.getPeak(getCurrentLayer().peaks).getIntensity();
			}
			else
			{
				PointType point = widgetToData_(last_mouse_pos_);
				mz = point[0] - start.getPeak(getCurrentLayer().peaks).getMZ();
				rt = point[1] - start.getSpectrum(getCurrentLayer().peaks).getRT();
				it = std::numeric_limits<DoubleReal>::quiet_NaN();
			}
		}
		else if (getCurrentLayer().type==LayerData::DT_CONSENSUS)
		{
			if (end.isValid())
			{
				mz = end.getFeature(getCurrentLayer().consensus).getMZ() - start.getFeature(getCurrentLayer().consensus).getMZ();
				rt = end.getFeature(getCurrentLayer().consensus).getRT() - start.getFeature(getCurrentLayer().consensus).getRT();
				it = end.getFeature(getCurrentLayer().consensus).getIntensity() / start.getFeature(getCurrentLayer().consensus).getIntensity();
			}
			else
			{
				PointType point = widgetToData_(last_mouse_pos_);
				mz = point[0] - start.getFeature(getCurrentLayer().consensus).getMZ();
				rt = point[1] - start.getFeature(getCurrentLayer().consensus).getRT();
				it = std::numeric_limits<DoubleReal>::quiet_NaN();
			}
		}
		else if (getCurrentLayer().type==LayerData::DT_CHROMATOGRAM)
		{
			//TODO CHROM
		}
		
		//draw text			
		QStringList lines;
		if (print_rt) lines.push_back("RT delta : " + QString::number(rt,'f',2));
		lines.push_back("m/z delta: " + QString::number(mz,'f',6));
		if (boost::math::isinf(it) || boost::math::isnan(it))
		{
			lines.push_back("int ratio: n/a");
		}
		else
		{
			lines.push_back("int ratio: " + QString::number(it,'f',2));			
		}
		drawText_(painter, lines);
	}

	void SpectrumCanvas::drawText_(QPainter& painter, QStringList text)
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
		for (int i=0;i<text.size(); ++i)
		{
			width = std::max(width, 4 + metrics.width(text[i]));
		}
		
		//draw backgrond for text
		painter.fillRect(2,3,width, height, QColor(255,255,255,200));
		
		//draw text
		painter.setPen(Qt::black);
		for (int i=0;i<text.size(); ++i)
		{
			painter.drawText(3, 3 + (i+1) * line_spacing, text[i]);
		}
		
		painter.restore();
	}

} //namespace

