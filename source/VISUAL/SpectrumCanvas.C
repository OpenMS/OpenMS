// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>

// custom cursor for translating
#include "ICONS/handopen.xpm"
#include "ICONS/handclosed.xpm"

// QT
#include <QtGui/QPainter>
#include <QtGui/QPaintEvent>
#include <QtGui/QBitmap>

using namespace std;

namespace OpenMS
{	
	SpectrumCanvas::SpectrumCanvas(QWidget* parent)
		: QWidget(parent),
			buffer_(),
			action_mode_(AM_ZOOM),
			intensity_mode_(IM_NONE),
			layers_(),
			mz_to_x_axis_(true),
			visible_area_(AreaType::empty),
			overall_data_range_(DRange<3>::empty),
			show_grid_(true),
			show_reduced_(false),
			update_buffer_(false),
			cursor_translate_(QPixmap(XPM_handopen)),
			cursor_translate_in_progress_(QPixmap(XPM_handclosed)),
			current_layer_(0),
			spectrum_widget_(0),
			datareducer_(0),
			percentage_factor_(1.0),
			snap_factor_(1.0),
			rubber_band_(QRubberBand::Rectangle,this)
	{
		setAttribute(Qt::WA_OpaquePaintEvent);
		// get mouse coordinates while mouse moves over diagramm.	
		setMouseTracking(TRUE);
		// prevents errors caused by too small width,height values
		setMinimumSize(200,200);
		// Take as much space as possible
		setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
	  
	  //reserve enough space to avoid copying layer data
	  layers_.reserve(10);
	}
	
	void SpectrumCanvas::resizeEvent(QResizeEvent* /* e */)
	{
#ifdef DEBUG_TOPPVIEW
		cout << "BEGIN " << __PRETTY_FUNCTION__ << endl;
#endif
		buffer_ = QPixmap(width(), height());
		update_buffer_ = true;
		updateScrollbars_();
		update_(__PRETTY_FUNCTION__);
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
	}
	
	void SpectrumCanvas::setDispInt(float min, float max)
	{
		layers_[current_layer_].min_int = min;
		layers_[current_layer_].max_int = max;
		intensityDistributionChange_();
	}
	
	void SpectrumCanvas::showGridLines(bool show)
	{
		show_grid_ = show;
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	void SpectrumCanvas::intensityDistributionChange_()
	{
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	void SpectrumCanvas::intensityModeChange_()
	{
		recalculateSnapFactor_();
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}

	bool SpectrumCanvas::isMzToXAxis() 
	{ 
		return mz_to_x_axis_; 
	}
	
	void SpectrumCanvas::mzToXAxis(bool mz_to_x_axis)
	{
		mz_to_x_axis_ = mz_to_x_axis;
		axisMappingChange_();
	}

	void SpectrumCanvas::axisMappingChange_()
	{
		updateScrollbars_();
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}

	void SpectrumCanvas::actionModeChange_()
	{
		
	}
	
	void SpectrumCanvas::changeVisibleArea_(const AreaType& new_area, bool add_to_stack)
	{
#ifdef DEBUG_TOPPVIEW
		cout << "BEGIN " << __PRETTY_FUNCTION__ << endl;
#endif
		if (new_area==visible_area_)
		{
			return;
		}
		//store old zoom state
		if (add_to_stack)
		{
			zoom_stack_.push(visible_area_);
		}
		visible_area_ = new_area;
		
		updateScrollbars_();
	
		emit visibleAreaChanged(new_area);
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
	}
	
	void SpectrumCanvas::updateScrollbars_()
	{
		
	}
	
	void SpectrumCanvas::zoomBack_()
	{
		if (zoom_stack_.empty())
		{
			resetZoom();
		}
		else
		{
			changeVisibleArea_(zoom_stack_.top());
			zoom_stack_.pop();
		}
	}
	
	void SpectrumCanvas::resetZoom()
	{
#ifdef DEBUG_TOPPVIEW
		cout << "BEGIN " << __PRETTY_FUNCTION__ << endl;
#endif
		zoom_stack_ = stack<AreaType>();

		AreaType tmp;
		tmp.assign(overall_data_range_);
		changeVisibleArea_(tmp);
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
	}
	
	void SpectrumCanvas::setVisibleArea(AreaType area)
	{
		changeVisibleArea_(area);
	}
	
	
	void SpectrumCanvas::paintGridLines_(QPainter& painter)
	{
		QColor g1(130,130,130);
		QColor g2(170,170,170);
		QColor g3(230,230,230);
		
		if (!show_grid_ || !spectrum_widget_) return;
	
		painter.save();
	
		unsigned int xl, xh, yl, yh; //width/height of the diagram area, x, y coordinates of lo/hi x,y values
	
		xl = 0;
		xh = width();

		yl = height();
		yh = 0;

	
		// drawing of grid lines and associated text	
		for (unsigned int j = 0; j != spectrum_widget_->xAxis()->gridLines().size() ; j++) 
		{
			// style definitions
			switch(j)
			{
				case 0:	// style settings for big intervals 
					painter.setPen(QPen(g1));
					break;
				case 1:	// style settings for small intervals
					painter.setPen(QPen(g2));
					break;
				case 2: // style settings for smalles intervals
					painter.setPen(QPen(g3));
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
		
		for (unsigned int j = 0; j != spectrum_widget_->yAxis()->gridLines().size() ; j++) 
		{

			// style definitions
			switch(j)
			{
				case 0:	// style settings for big intervals 
					painter.setPen(QPen(g1));
					break;
				case 1:	// style settings for small intervals
					painter.setPen(QPen(g2));
					break;
				case 2: // style settings for smalles intervals
					painter.setPen(QPen(g3));
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
	
	void SpectrumCanvas::setMainPreferences(const Param& prefs)
	{
		prefs_ = prefs;
	}
	
	UnsignedInt SpectrumCanvas::activeLayerIndex() const
	{
		return current_layer_;	
	}

	SpectrumCanvas::ExperimentType& SpectrumCanvas::addEmptyPeakLayer()
	{
		UnsignedInt newcount = getLayerCount()+1;
		layers_.resize(newcount);
		layers_.back().type = LayerData::DT_PEAK;
		return layers_[newcount-1].peaks;
	}

	SignedInt SpectrumCanvas::addLayer(const ExperimentType& in)
	{	
		layers_.resize(getLayerCount()+1);
		layers_.back().peaks = in;
		layers_.back().type = LayerData::DT_PEAK;
		return finishAdding();
	}

	SignedInt SpectrumCanvas::addLayer(const FeatureMapType& map, bool pairs)
	{
		layers_.resize(layers_.size()+1);
		layers_.back().features = map;
		if (pairs)
		{
			layers_.back().type = LayerData::DT_FEATURE_PAIR;
		}
		else
		{
			layers_.back().type = LayerData::DT_FEATURE;
		}
		return finishAdding();
	}
	
	void SpectrumCanvas::changeVisibility(int i , bool b)
	{
		LayerData& layer = getLayer_(i);
		if (layer.visible!=b)
		{
			layer.visible=b;
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
	}

  const DRange<3>& SpectrumCanvas::getDataRange()
  {
  	return overall_data_range_;
  }

	void SpectrumCanvas::updateRanges_(UnsignedInt layer_index, UnsignedInt mz_dim, UnsignedInt rt_dim, UnsignedInt it_dim)
	{
		if (layer_index >= getLayerCount())
		{
			return;
		}
		
		//update mz and RT
		DRange<3>::PositionType min = overall_data_range_.min();
		DRange<3>::PositionType max = overall_data_range_.max();
		
		if (getLayer(layer_index).type==LayerData::DT_PEAK)
		{
			if(show_reduced_)
			{
				const ExperimentType& red = getLayer(layer_index).reduced;
				if (red.getMinMZ() < min[mz_dim]) min[mz_dim] = red.getMinMZ();
				if (red.getMaxMZ() > max[mz_dim]) max[mz_dim] = red.getMaxMZ();
				if (red.getMinRT() < min[rt_dim]) min[rt_dim] = red.getMinRT();
				if (red.getMaxRT() > max[rt_dim]) max[rt_dim] = red.getMaxRT();
				if (red.getMinInt() < min[it_dim]) min[it_dim] = red.getMinInt();
				if (red.getMaxInt() > max[it_dim]) max[it_dim] = red.getMaxInt();
			}
			else
			{
				const ExperimentType& peaks = getLayer(layer_index).peaks;
				if (peaks.getMinMZ() < min[mz_dim]) min[mz_dim] = peaks.getMinMZ();
				if (peaks.getMaxMZ() > max[mz_dim]) max[mz_dim] = peaks.getMaxMZ();
				if (peaks.getMinRT() < min[rt_dim]) min[rt_dim] = peaks.getMinRT();
				if (peaks.getMaxRT() > max[rt_dim]) max[rt_dim] = peaks.getMaxRT();
				if (peaks.getMinInt() < min[it_dim]) min[it_dim] = peaks.getMinInt();
				if (peaks.getMaxInt() > max[it_dim]) max[it_dim] = peaks.getMaxInt();
		
			}
		}
		else
		{
			const FeatureMapType& feat = getLayer(layer_index).features;
			if (feat.getMin()[1] < min[mz_dim]) min[mz_dim] = feat.getMin()[1];
			if (feat.getMax()[1] > max[mz_dim]) max[mz_dim] = feat.getMax()[1];
			if (feat.getMin()[0] < min[rt_dim]) min[rt_dim] = feat.getMin()[0];
			if (feat.getMax()[0] > max[rt_dim]) max[rt_dim] = feat.getMax()[0];
			if (feat.getMinInt() < min[it_dim]) min[it_dim] = feat.getMinInt();
			if (feat.getMaxInt() > max[it_dim]) max[it_dim] = feat.getMaxInt();
		}
		
		overall_data_range_.setMin(min);
		overall_data_range_.setMax(max);
		
		//cout << "Updated range: " << overall_data_range_ << endl;
	}
	
	void SpectrumCanvas::recalculateRanges_(UnsignedInt mz_dim, UnsignedInt rt_dim, UnsignedInt it_dim)
	{
		overall_data_range_ = DRange<3>::empty;
		
		for (UnsignedInt i=0; i< getLayerCount(); ++i)
		{
			updateRanges_(i, mz_dim, rt_dim, it_dim);
		}
	}

	double SpectrumCanvas::getSnapFactor()
	{
		return snap_factor_;
	}

	void SpectrumCanvas::repaintAll()
	{
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
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
	
} //namespace


