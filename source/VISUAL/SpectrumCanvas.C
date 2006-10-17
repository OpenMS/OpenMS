// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
#include <qpainter.h>
#include <qpixmap.h>
#include <qbitmap.h>

// STLasdsd
#include <iostream>


using namespace std;
namespace OpenMS
{

	SpectrumCanvas::SpectrumCanvas(QWidget* parent, const char* name, WFlags f)
		: QWidget(parent, name, f | WNoAutoErase | WStaticContents ),
			buffer_(0),
			action_mode_(AM_ZOOM),
			intensity_mode_(IM_NONE),
			disp_ints_(),
			mz_to_x_axis_(true),
			pen_width_(0),
			show_grid_(true),
			recalculate_(false),
			current_data_(0),
			spectrum_widget_(0),
			datasets_(),
			features_(),
			percentage_factor_(1.0),
			snap_factor_(1.0)
	{
		// get mouse coordinates while mouse moves over diagramm.	
		setMouseTracking(TRUE);
		// prevents errors caused by too small width,height values
		setMinimumSize(200,200);
		// Take as much space as possible
		setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
		resetRanges_();
	  createCustomMouseCursors_();
	  
	  //reserve enough space for 20 datasets
	  datasets_.reserve(20);
	  features_.reserve(20);
	  
		// we need to initialize the painting buffer in the
		// constructor for maximum performance
		adjustBuffer_();
	}
	
	void SpectrumCanvas::createCustomMouseCursors_()
	{
		// create custom mouse cursor for translate action as Qt doesn't provide one
		QPixmap* pm1 = new QPixmap(handopen);
		pm1->setMask(pm1->createHeuristicMask());
		cursor_translate_ = QCursor(*pm1);
	
		QPixmap* pm2 = new QPixmap(handclosed);
		pm2->setMask(pm2->createHeuristicMask());
		cursor_translate_in_progress_ = QCursor(*pm2); 
	}
	
	void SpectrumCanvas::resizeEvent(QResizeEvent* /* e */)
	{
		adjustBuffer_();
		updateScrollbars_();
		invalidate_();
	}
	
	void SpectrumCanvas::paintEvent(QPaintEvent* e)
	{
		// blit changed parts of backbuffer on widget
		QMemArray<QRect> rects = e->region().rects();
		for (int i = 0; i < (int)rects.size(); ++i)
		{
			bitBlt(this, rects[i].topLeft(), buffer_, rects[i]);
		}
	}
	
	void SpectrumCanvas::adjustBuffer_()
	{
		if (buffer_) delete buffer_;
		
		buffer_ = new QPixmap(width(), height());
	}
	
	void SpectrumCanvas::setDispInt(float min, float max)
	{
		disp_ints_[current_data_] = make_pair(min,max);
		intensityDistributionChange_();
	}
	
	void SpectrumCanvas::showGridLines(bool show)
	{
		show_grid_ = show;
		invalidate_();
	}
	
	void SpectrumCanvas::intensityDistributionChange_()
	{
		recalculate_ = true;
		invalidate_();
	}
	
	void SpectrumCanvas::intensityModeChange_()
	{
		recalculateSnapFactor_();
		recalculate_ = true;
		invalidate_();
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
		recalculate_ = true;
		invalidate_();
	}

	void SpectrumCanvas::actionModeChange_()
	{
		
	}
	
	void SpectrumCanvas::changeVisibleArea_(const AreaType& new_area, bool add_to_stack)
	{
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
		recalculate_ = true;
		invalidate_();
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
		zoom_stack_ = stack<AreaType>();

		AreaType tmp;
		tmp.assign(overall_data_range_);
		changeVisibleArea_(tmp);
	}
	
	void SpectrumCanvas::setVisibleArea(AreaType area)
	{
		changeVisibleArea_(area);
	}
	
	
	void SpectrumCanvas::paintGridLines_(QPainter* p)
	{
		QColor g1(130,130,130);
		QColor g2(170,170,170);
		QColor g3(230,230,230);
		
		if (!show_grid_ || !spectrum_widget_) return;
	
		p->save();
	
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
					p->setPen(QPen(g1,pen_width_));
					break;
				case 1:	// style settings for small intervals
					p->setPen(QPen(g2,pen_width_));
					break;
				case 2: // style settings for smalles intervals
					p->setPen(QPen(g3,pen_width_));
					break;
				default:
					std::cout << "empty vertical grid line vector error!" << std::endl;
					p->setPen(QPen(QColor(0,0,0),pen_width_));
					break;
			}

			int x;
			for (std::vector<double>::const_iterator it = spectrum_widget_->xAxis()->gridLines()[j].begin(); it != spectrum_widget_->xAxis()->gridLines()[j].end(); it++) 
			{
				x = static_cast<int>(intervalTransformation(*it, spectrum_widget_->xAxis()->getAxisMinimum(), spectrum_widget_->xAxis()->getAxisMaximum(), xl, xh));
				p->drawLine(x, yl, x, yh);
			}
		}
		
		for (unsigned int j = 0; j != spectrum_widget_->yAxis()->gridLines().size() ; j++) 
		{

			// style definitions
			switch(j)
			{
				case 0:	// style settings for big intervals 
					p->setPen(QPen(g1,pen_width_));
					break;
				case 1:	// style settings for small intervals
					p->setPen(QPen(g2,pen_width_));
					break;
				case 2: // style settings for smalles intervals
					p->setPen(QPen(g3,pen_width_));
					break;
				default:
					std::cout << "empty vertical grid line vector error!" << std::endl;
					p->setPen(QPen(QColor(0,0,0),pen_width_));
					break;
			}

			int y;
			for (std::vector<double>::const_iterator it = spectrum_widget_->yAxis()->gridLines()[j].begin(); it != spectrum_widget_->yAxis()->gridLines()[j].end(); it++) 
			{
				y = static_cast<int>(intervalTransformation(*it, spectrum_widget_->yAxis()->getAxisMinimum(), spectrum_widget_->yAxis()->getAxisMaximum(), yl, yh));
				p->drawLine(xl, y, xh, y);
			}
		}
		
		p->restore();
	}
	
	void SpectrumCanvas::setMainPreferences(const Param& prefs)
	{
		prefs_ = prefs;
	}
	
	UnsignedInt SpectrumCanvas::activeDataSetIndex() const
	{
		return current_data_;	
	}
	
	UnsignedInt SpectrumCanvas::getDataSetCount() const
	{
 		 if(show_reduced_)
			{
  			return reduced_datasets_.size();
			}
		 
		return datasets_.size();
	}
		
	const String& SpectrumCanvas::getDataSetName(UnsignedInt index) const
	{
 		if(show_reduced_)
 		{
			reduced_datasets_[index].getName();
 		}
		
		if (index >= getDataSetCount())
		{
			return String::EMPTY;
		}
		if (type_[index] == DT_PEAK)
		{
			return datasets_[index].getName();
		}
		else
		{
			return features_[index].getName();
		}
	}
	
	bool SpectrumCanvas::isDataSetVisible(UnsignedInt index) const
	{
		return layer_visible_[index];
	}		

	SpectrumCanvas::ExperimentType& SpectrumCanvas::addEmptyDataSet()
	{
		datasets_.push_back(ExperimentType());
		features_.push_back(FeatureMapType());
		type_.push_back(DT_PEAK);
		return datasets_[getDataSetCount()-1];
	}

	SpectrumCanvas::ExperimentType& SpectrumCanvas::getDataSet_(UnsignedInt index) throw (Exception::IndexOverflow)
	{
		if (index >= getDataSetCount())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,index,getDataSetCount()-1);
		}
		if(show_reduced_)
		{
			return reduced_datasets_[index];
		}
		return datasets_[index];
	}

	SpectrumCanvas::ExperimentType& SpectrumCanvas::currentDataSet_() throw (Exception::IndexOverflow)
	{
		if (getDataSetCount()==0)
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,0,0);
		}
		if(show_reduced_)
		{
			return reduced_datasets_[current_data_];
		}
		return datasets_[current_data_];
	}

	const SpectrumCanvas::FeatureMapType& SpectrumCanvas::currentFeatureMap() const throw (Exception::IndexOverflow)
	{
		if (getDataSetCount()==0)
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,0,0);
		}
		return features_[current_data_];
	}

	SpectrumCanvas::FeatureMapType& SpectrumCanvas::currentFeatureMap_() throw (Exception::IndexOverflow)
	{
		if (getDataSetCount()==0)
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,0,0);
		}
		return features_[current_data_];
	}

	const SpectrumCanvas::ExperimentType& SpectrumCanvas::getDataSet(UnsignedInt index) const throw (Exception::IndexOverflow)
	{
		if (index >= getDataSetCount())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,index,getDataSetCount()-1);
		}
		if(show_reduced_)
		{
			return reduced_datasets_[index];
		}
		return datasets_[index];
	}

	const SpectrumCanvas::ExperimentType& SpectrumCanvas::currentDataSet() const throw (Exception::IndexOverflow)
	{
		if (getDataSetCount()==0)
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,0,0);
		}
		if(show_reduced_)
		{
			return reduced_datasets_[current_data_];
		}
		return datasets_[current_data_];
	}

	SignedInt SpectrumCanvas::addDataSet(const ExperimentType& in)
	{	
		datasets_.push_back(in);
		features_.push_back(FeatureMapType());
		type_.push_back(DT_PEAK);
		return finishAdding();
	}

	SignedInt SpectrumCanvas::addDataSet(const FeatureMapType& in)
	{
		datasets_.push_back(ExperimentType());	
		features_.push_back(in);
		type_.push_back(DT_FEATURE);
		return finishAdding();
	}
	
	void SpectrumCanvas::changeVisibility(int i , bool b)
	{
		if (layer_visible_[i]!=b)
		{
			layer_visible_[i]=b;
			invalidate_();
		}
	}

  const DRange<3>& SpectrumCanvas::getDataRange()
  {
  	return overall_data_range_;
  }

	void SpectrumCanvas::updateRanges_(UnsignedInt data_set, UnsignedInt mz_dim, UnsignedInt rt_dim, UnsignedInt it_dim)
	{
		
		if (data_set >= getDataSetCount())
		{
			return;
		}
		
		
		//update mz and RT
		DRange<3>::PositionType min = overall_data_range_.min();
		DRange<3>::PositionType max = overall_data_range_.max();
		if(show_reduced_)
		{
			if (type_[data_set]==DT_PEAK)
					{
						if (reduced_datasets_[data_set].getMinMZ() < min[mz_dim]) min[mz_dim] = reduced_datasets_[data_set].getMinMZ();
						if (reduced_datasets_[data_set].getMaxMZ() > max[mz_dim]) max[mz_dim] = reduced_datasets_[data_set].getMaxMZ();
						if (reduced_datasets_[data_set].getMinRT() < min[rt_dim]) min[rt_dim] = reduced_datasets_[data_set].getMinRT();
						if (reduced_datasets_[data_set].getMaxRT() > max[rt_dim]) max[rt_dim] = reduced_datasets_[data_set].getMaxRT();
						if (reduced_datasets_[data_set].getMinInt() < min[it_dim]) min[it_dim] = reduced_datasets_[data_set].getMinInt();
						if (reduced_datasets_[data_set].getMaxInt() > max[it_dim]) max[it_dim] = reduced_datasets_[data_set].getMaxInt();
					}
			
		}
		else
			{
				if (type_[data_set]==DT_PEAK)
					{
						if (datasets_[data_set].getMinMZ() < min[mz_dim]) min[mz_dim] = datasets_[data_set].getMinMZ();
						if (datasets_[data_set].getMaxMZ() > max[mz_dim]) max[mz_dim] = datasets_[data_set].getMaxMZ();
						if (datasets_[data_set].getMinRT() < min[rt_dim]) min[rt_dim] = datasets_[data_set].getMinRT();
						if (datasets_[data_set].getMaxRT() > max[rt_dim]) max[rt_dim] = datasets_[data_set].getMaxRT();
						if (datasets_[data_set].getMinInt() < min[it_dim]) min[it_dim] = datasets_[data_set].getMinInt();
						if (datasets_[data_set].getMaxInt() > max[it_dim]) max[it_dim] = datasets_[data_set].getMaxInt();
					}
				else
					{
						if (features_[data_set].getMin()[1] < min[mz_dim]) min[mz_dim] = features_[data_set].getMin()[1];
						if (features_[data_set].getMax()[1] > max[mz_dim]) max[mz_dim] = features_[data_set].getMax()[1];
						if (features_[data_set].getMin()[0] < min[rt_dim]) min[rt_dim] = features_[data_set].getMin()[0];
						if (features_[data_set].getMax()[0] > max[rt_dim]) max[rt_dim] = features_[data_set].getMax()[0];
						if (features_[data_set].getMinInt() < min[it_dim]) min[it_dim] = features_[data_set].getMinInt();
						if (features_[data_set].getMaxInt() > max[it_dim]) max[it_dim] = features_[data_set].getMaxInt();
					}
			}
		overall_data_range_.setMin(min);
		overall_data_range_.setMax(max);
		
		//cout << "Updated range: " << overall_data_range_ << endl;
	}

	void SpectrumCanvas::resetRanges_()
	{
		overall_data_range_ = DRange<3>::empty;
		//cout << "Reset range: " << overall_data_range_ << endl;
	}
	
	void SpectrumCanvas::recalculateRanges_(UnsignedInt mz_dim, UnsignedInt rt_dim, UnsignedInt it_dim)
	{
		resetRanges_();
		
		for (UnsignedInt i=0; i< getDataSetCount(); ++i)
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
		recalculate_ = true;
		invalidate_();
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
	
	SpectrumCanvas::DataType SpectrumCanvas::getCurrentDataType() const
	{
		return type_[current_data_];
	}
		
} //namespace


