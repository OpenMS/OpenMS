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
// $Id: Spectrum1DCanvas.C,v 1.53 2006/06/09 22:00:08 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// Qt
#include <qmessagebox.h>
#include <qaction.h>
#include <qtimer.h>

// OpenMS
#include <OpenMS/VISUAL/PeakIcon.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DCanvasPDP.h>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	using namespace Internal;
	
	Spectrum1DCanvas::Spectrum1DCanvas(QWidget* parent, const char* name, WFlags f)
		: SpectrumCanvas(parent, name, f | WRepaintNoErase),
		zoom_factor_(2.0),
		absolute_intensity_(false),
		layer_factor_(1.0),
		max_layer_(0),
		snap_to_max_mode_(false),
		snap_factor_(1.0),
		data_area_(1.0, 1.0, 1000.0, 100.0),
		zoom_status_(),
		is_highlighted_(false)
	{
		// set style definitions
		style_.border = 7;
	
		// get mouse coordinates while mouse moves over diagramm.	
		viewport()->setMouseTracking(TRUE);
	
		activateDataSet(0);
	}
	
	//change the current data set
	void Spectrum1DCanvas::activateDataSet(int data_set)
	{
		if (data_set >= int(getDataSetCount()))
		{
			return ;
		}
		
		current_data_ = data_set;
			
		// no peak is selected
		nearest_peak_ = currentDataSet()[0].end();
		selected_peaks_.clear();
		selected_peaks_.push_back(currentDataSet()[0].begin());
		
		emit layerActivated(this);
	}
	
	void Spectrum1DCanvas::viewportPaintEvent(QPaintEvent* e)
	{
		SpectrumCanvas::viewportPaintEvent(e);
		
		QPainter painter(viewport());
		rubber_band_.draw(painter);
		painter.end();
	}
	
	void Spectrum1DCanvas::setVisibleArea(double lo, double hi)
	{
		changeVisibleArea_(lo, hi);
	}
	
	void Spectrum1DCanvas::setVisibleArea(DRange<2> range)
	{
		changeVisibleArea_(range.minX(), range.maxX());
	}
	
	void Spectrum1DCanvas::changeVisibleArea_(double lo, double hi)
	{
		changeVisibleArea_(AreaType(lo, visible_area_.minY(), hi, visible_area_.maxY()));
	}
	
	QPoint Spectrum1DCanvas::chartToWidget_(const PeakType& peak)
	{
		return SpectrumCanvas::chartToWidget_(PointType(peak.getPosition()[0], snap_factor_*layer_factor_*peak.getIntensity()));
	}
	
	const Spectrum1DCanvas::AreaType& Spectrum1DCanvas::getDataArea_()
	{
		return data_area_;
	}
	
	void Spectrum1DCanvas::zoomIn(double position, int /*steps*/)
	{
		double delta, newLo, newHi;
	
		delta = (visible_area_.maxX()-visible_area_.minX())/(2*zoom_factor_);
		newLo = position - delta;
		newHi = position + delta;
	
	// 	setVisibleAreaAnimated(newLo, newHi, steps);
		changeVisibleArea_(newLo, newHi);
	}
	
	
	void Spectrum1DCanvas::zoomOut(double position,int /*steps*/)
	{
		double delta;
		double newLo;
		double newHi;
	
		delta = zoom_factor_*(visible_area_.maxX()-visible_area_.minX())/2;
		newLo = position-delta;
		newHi = position+delta;
	
		// fit too data set bounds
		if (newLo < data_area_.minX()) newLo = data_area_.minX();
		if (newHi > data_area_.maxX()) newHi = data_area_.maxX();
	
	// 	setVisibleAreaAnimated(newLo, newHi, steps);
		changeVisibleArea_(newLo, newHi);
	}
	
	void Spectrum1DCanvas::translate(double position,int /*steps*/)
	{
		double delta;
		double newLo;
		double newHi;
	
		delta = (visible_area_.maxX()-visible_area_.minX())/2;
		newLo = position-delta;
		newHi = position+delta;
	
		if (newLo < data_area_.minX())
		{
			changeVisibleArea_(data_area_.minX(), data_area_.minX() + 2 * delta);
			return;
		}
	
		if (newHi >data_area_.maxX())
		{
			changeVisibleArea_(data_area_.maxX() - 2 * delta, data_area_.maxX());
			return;
		}
	
		changeVisibleArea_(newLo, newHi);
	}
	
	//////////////////////////////////////////////////////////////////////////////////
	// Qt events
	
	void Spectrum1DCanvas::contentsMousePressEvent( QMouseEvent *e)
	{
		// get mouse position in widget coordinates
		QPoint p = contentsToViewport(e->pos());
		action_start_pos_ = p;
		action_current_pos_ = p;
	
		if (e->button() == LeftButton)
		{
			PointType pos = widgetToChart_(p);
			if (e->state() & QMouseEvent::ShiftButton)
			{
				if (action_mode_ == AM_ZOOM)
				{
					rubber_band_.show();
					rubber_band_.setTopLeft(p);
					rubber_band_.setBottomRight(p);
					rubber_band_.updateRegion(viewport());
					zoom_status_ = QString("(%1,%2) => ").arg(pos.X(),0,'f',2).arg(pos.Y(),0,'f',2);
				}
			}
			else
			{
				if (action_mode_ == AM_ZOOM)
				{
					rubber_band_.show();	
					rubber_band_.setTopLeft(chartToWidget_(PointType(pos.X(), data_area_.maxY())));
					rubber_band_.setBottomRight(chartToWidget_(PointType(pos.X(), data_area_.minY())));
	
					rubber_band_.updateRegion(viewport());
					zoom_status_ = QString("%1 => ").arg(pos.X(),0,'f',2);
				}
			}
			if (action_mode_ == AM_TRANSLATE) viewport()->setCursor(cursor_translate_in_progress_);
		}
	}
	
	
	void Spectrum1DCanvas::contentsMouseReleaseEvent(QMouseEvent *e)
	{
		switch (action_mode_)
		{
			case AM_SELECT:
			{
				if (e->button() == Qt::RightButton)
				{
					emit contextMenu(e->globalPos());
				}
				// Peak selection
				else if (nearest_peak_ != currentDataSet()[0].end())
				{
					DataValue tmp = nearest_peak_->getMetaValue(4);
					if (tmp.isEmpty() || (!tmp.isEmpty() && UnsignedInt(tmp) == Spectrum1DCanvas::IT_NOICON))
					{
						nearest_peak_->setMetaValue(4, SignedInt(Spectrum1DCanvas::IT_CIRCLE));
						selected_peaks_.push_back(nearest_peak_);
	
						std::ostringstream msg;
						msg << "Selected peak at position " << nearest_peak_->getPosition()[0]  << " (" << (selected_peaks_.size()-1);
						msg << " peaks selected altogether.)";
						emit sendStatusMessage(msg.str(), 5000);
					}
					else
					{
						nearest_peak_->setMetaValue(4,SignedInt(Spectrum1DCanvas::IT_NOICON));
						std::vector<SpectrumIteratorType>::iterator it_tmp =
						  std::find(selected_peaks_.begin(), selected_peaks_.end(), nearest_peak_);
	
						if(it_tmp != selected_peaks_.end())
						{
							selected_peaks_.erase(it_tmp);
	
							std::ostringstream msg;
							msg << "Deselected peak at position " << nearest_peak_->getPosition()[0]  << " (" << (selected_peaks_.size()-1);
							msg << " peaks selected altogether.)";
							emit sendStatusMessage(msg.str(), 5000);
						}
	
						//std::cout << "selected_peaks_.size(): " << selected_peaks_.size() << std::endl;
					}
					invalidate_();
				}
				break;
			}
			case AM_ZOOM:
			{
				// zoom-in-at-position or zoom-in-to-area
				if (e->button() == LeftButton)
				{
					zoom_status_ = "";
					rubber_band_.hide();
					rubber_band_.updateRegion(viewport());
					QRect rect = rubber_band_.getRect();
					rubber_band_.setRect(QRect(0,0,0,0));
	
					if (rect.width() > 4 && rect.height() > 4)   // free zoom
					{
						Spectrum1DCanvas::AreaType area(widgetToChart_(rect.topLeft()), widgetToChart_(rect.bottomRight()));
						//area.normalize();
						if (e->state() & QMouseEvent::ShiftButton)
						{
							//check if selected area is outside data area and correct errors
							if (area.minX() < data_area_.minX())
							{
								area.setMinX(data_area_.minX());
							}
							if (area.minY() < data_area_.minY())
							{
								area.setMinY(data_area_.minY());
							}
							if (area.maxX() > data_area_.maxX())
							{
								area.setMaxX(data_area_.maxX());
							}
							if (area.maxY() > data_area_.maxY())
							{
								area.setMaxY(data_area_.maxY());
							}
	
							changeVisibleArea_(area);
						}
						else
						{
							changeVisibleArea_(area.minX(), area.maxX());
						}
					}
					else                                        // position axis only zoom
					{
						zoomIn(widgetToChart_(rect.topLeft()).X(), 10);
					}
				}
	
				// zoom-back
				if (e->button() == Qt::MidButton)
				{
					zoomBack_();
				}
	
				if (e->button() == Qt::RightButton)
				{
					emit contextMenu(e->globalPos());
				}
				break;
			}
			case AM_TRANSLATE:
			{
				viewport()->setCursor(cursor_translate_);
				if (e->button() == Qt::RightButton)
				{
					emit contextMenu(e->globalPos());
				}
				break;
			}
			default:/* throw new std::runtime_error("undefined action mode: action_mode_="+action_mode_);*/
				break;
		}
	}
	
	void Spectrum1DCanvas::contentsMouseDoubleClickEvent( QMouseEvent *e)
	{
		if (e->button() == LeftButton && action_mode_ == AM_SELECT)
		{
			SpectrumIteratorType i = findPeakAtPosition(contentsToViewport(e->pos()));
			if (i != currentDataSet()[0].end())
			{
				i->metaRegistry().registerName("extended_label","","");
				if (i->metaValueExists("extended_label"))
				{
					QMessageBox::information( this, "Extended meta information",
					QString(i->getMetaValue("extended_label").toChar()) );
				}
			}
		}
	
		// mid-doubleclick shows the whole spectrum
		if (e->button() == Qt::MidButton)
		{
			resetZoom();
		}
	}
	
	void Spectrum1DCanvas::contentsMouseMoveEvent( QMouseEvent *e)
	{
		// mouse position relative to the diagram widget
		QPoint p = contentsToViewport(e->pos());
	
		switch (action_mode_)
		{
			case AM_SELECT:
			{ 
				emit sendCursorStatus();
	 			viewport()->setCursor(Qt::ArrowCursor);
				is_highlighted_ = true;
				nearest_peak_ = findPeakAtPosition(p);
				invalidate_();
				break;
			}
			case AM_ZOOM:
			{
				is_highlighted_ = true;
				viewport()->setCursor(Qt::CrossCursor);
				PointType pos = widgetToChart_(p);
	
				if (e->state() & LeftButton)
				{
					if (e->state() & QMouseEvent::ShiftButton) // free zoom
					{
						rubber_band_.updateRegion(viewport());
						rubber_band_.setBottomRight(p);
						rubber_band_.updateRegion(viewport());
					}
					else // zoom on position axis only
					{
						rubber_band_.updateRegion(viewport());
						rubber_band_.setBottomRight(chartToWidget_(PointType(pos.X(), data_area_.minY())));
						rubber_band_.updateRegion(viewport());
					}
				}
				if (e->state() & QMouseEvent::ShiftButton)
					emit sendCursorStatus( pos.X(), pos.Y());
				else
					emit sendCursorStatus( pos.X() );
				break;
			}
			case AM_TRANSLATE:
			{
				viewport()->setCursor(cursor_translate_);
				// Translation of visible area
				if(e->state() & QMouseEvent::LeftButton)
				{
					if (action_current_pos_!=action_start_pos_)
					{
						viewport()->setCursor(cursor_translate_in_progress_);
						double x1,x2;
						double newLo, newHi;
						double dx;
						x1 = widgetToChart_(action_current_pos_).X();
						x2 = widgetToChart_(p).X();
						dx = x1 - x2;  // translation in chart metric
						newLo = visible_area_.minX() + dx;
						newHi = visible_area_.maxX() + dx;
						// check if we are falling out of bounds
						if (newLo < data_area_.minX())
						{
							changeVisibleArea_(data_area_.minX(), data_area_.minX()+(visible_area_.maxX()-visible_area_.minX()));
						}
						else if (newHi > data_area_.maxX())
						{
							changeVisibleArea_(data_area_.maxX()-(visible_area_.maxX()-visible_area_.minX()), data_area_.maxX());
						}
						else
						{
							changeVisibleArea_(newLo, newHi);
						}
					}
					action_current_pos_=p;
				// End of: Translation
				}
				break;
			}
			default: /*throw new std::runtime_error("undefined action mode: action_mode_="+action_mode_);*/
				break;
		}
	}
	
	Spectrum1DCanvas::SpectrumIteratorType Spectrum1DCanvas::findPeakAtPosition(QPoint p)
	//  Input: position p on screen
	//  Algorithmn:
	//	1. step:
	//		transform interval [P.x, P.x+1.0) into diagramm metrics: [P.x, P.x+1.0) => [d.x1, d.x2)
	//	2. step:
	//		find the 2 iterators in the data set that points to the first peak that falls into this interval
	//		and the first, that falls off the interval
	//	Depending on how many peaks (if any) lie between this two iterators further decisions have to be made
	//  Case 0: no peak
	//	Case 1: one peak => we already found the nearest peak
	// 	Case 2: several peaks => as all peaks between these two iterators will be projected on the same line on the screen further
	//					decision must me made (based on P.y) (see code)
	// now a slightly modified version is used. the only difference is, that the source interval
	// has a length of 3 pixels (to simplify grabbing of a peak).
	
	{
		// get the interval (in diagramm metric) that will be projected on screen coordinate p.x() or p.y() (depending on orientation)
	// 	int coord;
		double interval_start, interval_end;
		PointType lt = widgetToChart_(p - QPoint(1, 1));
		PointType rb = widgetToChart_(p + QPoint(1, 1));
	/*	if (mapping_info_.isMzToXAxis())
		{*/
			interval_start = lt.X();
			interval_end = rb.X();
	/*	}
		else
		{
			interval_start = lt.y;
			interval_end = rb.y;
		}*/
	
		// debug code:
		//	std::cout << "Intervall start: " << interval_start << std::endl;
		//	std::cout << "Intervall end:   " << interval_end << std::endl;
	
		// get iterator on first peak with higher position than interval_start
		PeakType temp;
		temp.getPosition()[0] = interval_start;
		SpectrumIteratorType left_it = std::lower_bound(visible_begin_[current_data_], visible_end_[current_data_], temp, PeakType::PositionLess());
	
		// get iterator on first peak with higher position than interval_end
		temp.getPosition()[0] = interval_end;
		SpectrumIteratorType	right_it = std::lower_bound(visible_begin_[current_data_], visible_end_[current_data_], temp, PeakType::PositionLess());
	
		//debug code:
		//	std::cout << "left_it and *left_it: "  << left_it->getPosition()[0] << " " << left_it->getIntensity() << std::endl;
		//	std::cout << "right_it and *right_it: "  << right_it->getPosition()[0] << " " << right_it->getIntensity() << std::endl;
	
		if (left_it == right_it)
		{
		//	std::cout << "case 0: no peak" << std::endl;	// debug code
			return currentDataSet()[0].end();  // both are equal => no peak falls into this interval
		}
	
		if (left_it == right_it-1 )
		{
		//	std::cout << "case 1: one peak falls into the interval " << std::endl; // debug code
			return left_it;
		}
	
		// debug code
		//	std::cout << " case 3: several peaks in the same interval" << std::endl;
		//	std::cout << left_it->getIntensity() << " " << right_it->getIntensity() << std::endl;
	
		SpectrumIteratorType nearest_it = left_it;
	
	// 	double dest_interval_start = chartToWidget_(Point<float>(0, data_area_.minY())).y();
	// 	double dest_interval_end = chartToWidget_(Point<float>(0, data_area_.maxY())).y();
		double dest_interval_start, dest_interval_end;
		// select source interval start and end depending on diagram orientation
	
		if (mapping_info_.isMzToXAxis())
		{
			// select destination interval start and end depending on diagram orientation AND axis direction
			if (mapping_info_.isYAxisAsc())
			{
				dest_interval_start = viewport()->height();
				dest_interval_end = 0;
			}
			else
			{
				dest_interval_start = 0;
				dest_interval_end = viewport()->height();
			}
		}
		else
		{
			// select destination interval start and end depending on diagram orientation and axis direction
			if (mapping_info_.isXAxisAsc())
			{
				dest_interval_start = 0;
				dest_interval_end = viewport()->width();
			}
			else
			{
				dest_interval_start = viewport()->width();
				dest_interval_end = 0;
			}
		}
	
		int nearest_intensity = static_cast<int>(intervalTransformation(nearest_it->getIntensity(), visible_area_.minY(),
		                                                                 visible_area_.maxY(), dest_interval_start, dest_interval_end));
		int current_intensity;
	
		for (SpectrumIteratorType it = left_it; it != right_it; it++)
		{
			current_intensity = static_cast<int>(intervalTransformation(it->getIntensity(), visible_area_.minY(), visible_area_.maxY(),
			                                                             dest_interval_start, dest_interval_end));
			if (mapping_info_.isMzToXAxis())
			{
				if ( abs(current_intensity - p.y()) < abs(nearest_intensity - p.y()))
				{
					nearest_intensity = current_intensity;
					nearest_it = it;
				}
			}
			else
			{
				if ( abs(current_intensity - p.x()) < abs(nearest_intensity - p.x()))
				{
					nearest_intensity = current_intensity;
					nearest_it = it;
				}
			}
		}
	
		return nearest_it;
	}
	
	//////////////////////////////////////////////////////////////////////////////////
	// SLOTS
	
	void Spectrum1DCanvas::removeDataSet(int data_set)
	{
		if (data_set >= int(getDataSetCount()))
		{
			return;
		}
	
		//remove the data
		datasets_.erase(datasets_.begin()+data_set);
	
		//remove visibility setting
		layer_visible_.erase(layer_visible_.begin()+data_set);
	
		//renew values from visible_begin_ and visible_end_
		visible_begin_.clear();
		visible_end_.clear();
		for (UnsignedInt index=0; index < getDataSetCount(); ++index)
		{
			visible_begin_.push_back(getDataSet(index)[0].begin());
			visible_end_.push_back(getDataSet(index)[0].end());
		}
	
		//remove draw mode
		draw_modes_.erase(draw_modes_.begin()+data_set);
	
		//update current data set
		if (getDataSetCount() >= current_data_)
		{
			--current_data_;
		}
	
		if (datasets_.empty())
		{
			return;
		}
	
		//update visible area
		data_area_.setMinX(getDataSet(0)[0].getMin()[0]);
		data_area_.setMaxX(getDataSet(0)[0].getMax()[0]);
		data_area_.setMinY(1.0);  // minimal intensity always 1.0
		data_area_.setMaxY(getDataSet(0)[0].getMaxInt());
		max_layer_ = 0;
	
		for (UnsignedInt index = 1; index < getDataSetCount(); ++index)
		{
			if (data_area_.minX() > getDataSet(index)[0].getMin()[0])
			{
				data_area_.setMinX(getDataSet(index)[0].getMin()[0]);
			}
			if (data_area_.maxX() < getDataSet(index)[0].getMax()[0])
			{
				data_area_.setMaxX(getDataSet(index)[0].getMax()[0]);
			}
	
			if (data_area_.maxY() < getDataSet(index)[0].getMaxInt())
			{
				data_area_.setMaxY(getDataSet(index)[0].getMaxInt());
				max_layer_ = index;
			}
		}
	
		// extend region for cosmectical reasons -> peak with smallest/highest position won't fall on diagramm border on maximized view
		data_area_.setMaxY(data_area_.maxY() + 0.002 * (data_area_.maxY() - data_area_.minY()));
		data_area_.setMinX(data_area_.minX() - 0.002 * (data_area_.maxX() - data_area_.minX()));
		data_area_.setMaxX(data_area_.maxX() + 0.002 * (data_area_.maxX() - data_area_.minX()));
	
	
		changeVisibleArea_(data_area_);
	
		//
		if (data_area_.maxX() - data_area_.minX() <1.0)
		{
			changeVisibleArea_(data_area_.minX() -1.0, data_area_.maxX() + 1.0);
		}
		else
		{
			changeVisibleArea_(data_area_.minX(), data_area_.maxX());
		}
		showGridLines(show_grid_);
	}
	
	void Spectrum1DCanvas::setZoomFactor(double d)
	{
		zoom_factor_ = d;
	}
	
	void Spectrum1DCanvas::setDrawMode(QAction* a)
	{
		QString name = a->name();
	
		if (name == "setPeakMode")
			draw_modes_[current_data_] = DM_PEAKS;
		else if (name == "setConnectedLineMode")
			draw_modes_[current_data_]= DM_CONNECTEDLINES;
	
		//throw new std::runtime_error("unknown QAction");
		invalidate_();
	}
	
	void Spectrum1DCanvas::drawModePeaks()
	{
		draw_modes_[current_data_] = DM_PEAKS;
		invalidate_();
	}
	
	void Spectrum1DCanvas::drawModeLines()
	{
		draw_modes_[current_data_] = DM_CONNECTEDLINES;
		invalidate_();
	}
	
	//////////////////////////////////////////////////////////////////////////////////
	// data plotting
	void Spectrum1DCanvas::drawPoints_(UnsignedInt index)
	{
		painter_.save();
		painter_.setBrush(NoBrush);
	
		QPoint p;
		for (SpectrumIteratorType i = visible_begin_[index]; i != visible_end_[index]; ++i)
		{
			// project diagramm coordinates on paint device coordinates
			p = chartToWidget_(*i);
	
			// set normal/highlighted pen
			painter_.setPen(i == nearest_peak_ ? high_pen_ : norm_pen_);
			// draw point
			painter_.drawLine(p.x(), p.y() - 1, p.x(), p.y() + 1);
			painter_.drawLine(p.x() - 1, p.y(), p.x() + 1, p.y());
		}
	
		painter_.restore();
	}
	
	
	void Spectrum1DCanvas::drawIcon(const PeakType& peak, const QPoint& p)
	{
		// draw icon associated to current peak
		if (peak.metaValueExists(4))
		{
			DataValue tmp = peak.getMetaValue(4);
			painter_.setPen(icon_pen_);
	
			UnsignedInt rect_half_size = 5;
			QRect rect(p.x() - rect_half_size, p.y() - rect_half_size, 2 * rect_half_size, 2 * rect_half_size);
	
			switch(UnsignedInt(tmp))
			{
				case Spectrum1DCanvas::IT_CIRCLE:
					PeakIcon::drawEllipse(painter_, rect);
					break;
				case Spectrum1DCanvas::IT_TRIANGLE:
					PeakIcon::drawTriangle(painter_, rect);
					break;
				case Spectrum1DCanvas::IT_ASTERIX:
					PeakIcon::drawAsterix(painter_, rect);
					break;
				case Spectrum1DCanvas::IT_SQUARE:
					PeakIcon::drawRectangle(painter_, rect);
					break;
				default:
					break;
			}
			painter_.setPen(norm_pen_);
		}
	}
	
	void Spectrum1DCanvas::drawPeaks_(UnsignedInt index)
	{
		// font for drawing of values on highlighted peaks
		QFont high_font(QFont("courier", static_cast<UnsignedInt>(1.5 * grid_row_width_)));
	
		painter_.save();
		painter_.setBrush(NoBrush);
	
		QPoint p, p0;
		for (SpectrumIteratorType i = visible_begin_[index]; i != visible_end_[index]; ++i)
		{
			if (i->getIntensity() < min_disp_ints_[current_data_] || i->getIntensity() > max_disp_ints_[current_data_])
			{
				continue;
			}
	
			p = chartToWidget_(*i);
			p0 = chartToWidget_(PointType(i->getPosition()[0], 0.0f));
	
			// highlight selected peak
			if (i == nearest_peak_)
			{
				painter_.setFont(high_font);
				painter_.setPen(high_pen_);
	
				painter_.drawLine(p0, p);
	
				// convert peak value to original value if intensity in logarithmic scale
				bool is_intensity_axis_percent = !isAbsoluteIntensity();
				double tmp_intens = (intensity_modification_ == IM_NONE) ?
					i->getIntensity() : log2linear(i->getIntensity(),is_intensity_axis_percent,old_max_intensity_);
	
				emit sendCursorStatus( i->getPosition()[0], tmp_intens);
			}
			else
			{
				// custom peak color
				if (i->metaValueExists(5))
				{
					QColor color = QColor(string(i->getMetaValue(5)).c_str());
					QPen pen(color, pen_width_);
					painter_.setPen(pen);
				}
	
				else
				{
					painter_.setPen(norm_pen_);
				}
				// draw normal peak
				painter_.drawLine(p0, p);
			}
	
			drawIcon(*i, p);
		}
		painter_.restore();
	}
	
	void Spectrum1DCanvas::drawConnectedLines_(UnsignedInt index)
	{
		painter_.save();
		painter_.setBrush(NoBrush);
		// drawing of the peaks
		painter_.setPen(norm_pen_);
		bool firstPoint=true;
		if (visible_begin_[index]==visible_end_[index])
		{
			// check cases where no peak at all is visible
			if (visible_begin_[index] == getDataSet(index)[0].end()) return;
			if (visible_end_[index] == getDataSet(index)[0].begin()) return;
	
			// draw line (clipping performed by Qt on both sides)
			painter_.drawLine(chartToWidget_(*(visible_begin_[index] - 1)),
			                  chartToWidget_(*visible_begin_[index]));
			return;
		}
	
		// connect peaks in visible area; (no clipping needed)
		QPoint p;
		for (SpectrumIteratorType it = visible_begin_[index]; it != visible_end_[index]; it++)
		{
			// connect lines
			if (firstPoint)
			{
				p = chartToWidget_(*it);
				painter_.moveTo(p);
				firstPoint = false;
			} else
			{
				p = chartToWidget_(*it);
				// custom peak color
				if (it->metaValueExists(5))
				{
					QColor color = QColor(string(it->getMetaValue(5)).c_str());
					QPen pen(color, pen_width_);
					painter_.setPen(pen);
				}
				painter_.lineTo(p);
			};
	
			// highlight selected peak
			if (it == nearest_peak_)
			{
				painter_.save();
				painter_.setPen(high_pen_);
				painter_.drawLine(p.x(), p.y()-4, p.x(), p.y()+4);
				painter_.drawLine(p.x()-4, p.y(), p.x()+4, p.y());
				painter_.restore();
				double tmp_intens;
				if (intensity_modification_ == IM_NONE)
					tmp_intens = it->getIntensity();
				else
					tmp_intens = log2linear(it->getIntensity(),!isAbsoluteIntensity(),old_max_intensity_);
	
				emit sendCursorStatus( it->getPosition()[0], tmp_intens);
			}
			// draw associated icon
			drawIcon(*it, p);
		}
	
		// clipping on left side
		if (visible_begin_[index] > getDataSet(index)[0].begin())
		{
			painter_.drawLine(chartToWidget_(*(visible_begin_[index]-1)), chartToWidget_(*(visible_begin_[index])));
		}
	
		// clipping on right side
		if (visible_end_[index] < getDataSet(index)[0].end())
		{
			painter_.drawLine( chartToWidget_(*(visible_end_[index]-1)), chartToWidget_(*(visible_end_[index])));
		}
	
		painter_.restore();
	}
	
	void Spectrum1DCanvas::setMainPreferences(const Param& prefs)
	{
		SpectrumCanvas::setMainPreferences(prefs);
		mapping_info_.setParam(prefs.copy("Preferences:1D:Mapping:",true));
	}
	
	void Spectrum1DCanvas::setBounds_()
	{
		if (!datasets_.empty())
		{
			// get iterators on peaks that outline the visible area
			PeakType temp;
			temp.getPosition()[0] = visible_area_.minX();
			for (UnsignedInt i=0; i<getDataSetCount();++i)
			{
				visible_begin_[i] = std::upper_bound(getDataSet(i)[0].begin(), getDataSet(i)[0].end(), temp, PeakType::PositionLess());
			}
	
			temp.getPosition()[0] = visible_area_.maxX();
			for (UnsignedInt i=0; i<getDataSetCount();++i)
			{
				visible_end_[i] = std::upper_bound(visible_begin_[i], getDataSet(i)[0].end(), temp, PeakType::PositionLess());
			}
	
			// If snap-to-max-mode is on: find local max and set factor to increase all data appropriately
			if (snap_to_max_mode_) {
				double local_max  = std::max_element(visible_begin_[0], visible_end_[0], PeakType::IntensityLess())->getIntensity();
				for (UnsignedInt i=1; i<getDataSetCount();++i)
				{
					SpectrumIteratorType tmp  = std::max_element(visible_begin_[i], visible_end_[i], PeakType::IntensityLess());
					if (tmp->getIntensity() > local_max) local_max = tmp->getIntensity();
				}
				snap_factor_ = 1.0*datasets_[max_layer_][0].getMaxInt()/local_max;
			}
			else 
				snap_factor_ = 1.0;
	
			if (action_mode_ != AM_SELECT)
				nearest_peak_ = visible_end_[current_data_];
		}
	}
	
	void Spectrum1DCanvas::invalidate_()
	{
		//cout << "INVALIDATE"<<endl;
		setBounds_();
	
		// get color settings
		setPaletteBackgroundColor(QColor(getPrefAsString("Preferences:1D:BackgroundColor").c_str()));
		viewport()->setBackgroundColor(QColor(getPrefAsString("Preferences:1D:BackgroundColor").c_str()));
	
		norm_pen_ = QPen(QColor(getPrefAsString("Preferences:1D:PeakColor").c_str()), pen_width_);
		high_pen_ = QPen(QColor(getPrefAsString("Preferences:1D:HighColor").c_str()), pen_width_+2);
		icon_pen_ = QPen(QColor(getPrefAsString("Preferences:1D:IconColor").c_str()), pen_width_);
	
	// repaint content
		painter_.begin(buffer_);
		// flushing buffer_ with widget color
		painter_.fillRect(0, 0, viewport()->width(), viewport()->height(), viewport()->backgroundColor());
	
		paintGridLines_(&painter_);
	
		for (UnsignedInt i=0; i< getDataSetCount();++i)
		{
			if (layer_visible_[i]==false)
			{
				continue;
			}
	
			// If multiple layers are shown and drawn in relative scale: 
			// increase each spectrum maximum to 100%
			// TODO: relative log-scale
			if (isAbsoluteIntensity())
				layer_factor_ = 1.0;
			else 
				layer_factor_ = 1.0*datasets_[max_layer_][0].getMaxInt()/getDataSet(i)[0].getMaxInt();
	
			switch (draw_modes_[i])
			{
				case DM_PEAKS:
					drawPeaks_(i);
				break;
				case DM_CONNECTEDLINES:
					drawConnectedLines_(i);
				break;
				default: /*throw new std::runtime_error("undefined draw mode: draw_modes_="+draw_modes_);*/ break;
			}
		}
	
		painter_.end();
		viewport()->repaint(false);
		//cout << "/INVALIDATE"<<endl;
	}
	
	
	void Spectrum1DCanvas::intensityModificationChange_()
	{
		//cout << "IM_CHANGE" <<endl;
		if (intensity_modification_ == IM_LOG)
			scaleAllData_(true);
		else if (intensity_modification_ == IM_NONE)
			scaleAllData_(false);
		//cout << "/IM_CHANGE" <<endl;
		invalidate_();
	}
	
	void Spectrum1DCanvas::legendModificationChange_()
	{
		update();
	}
	
	void Spectrum1DCanvas::scaleData_(bool is_log)
	{
		//abort if there is no data
		if (getDataSetCount()==0 || currentDataSet()[0].size()==0)	return ;
	
		bool is_intensity_axis_percent = !isAbsoluteIntensity();
	
		if (is_log)
		{
			old_max_intensity_ = data_area_.maxY();
			for (SpectrumIteratorType it=currentDataSet()[0].begin();it!=currentDataSet()[0].end();++it)
			{
				if (it->getIntensity()!=0)
				{
					//				cerr<<  it->getIntensity() << " => ";
					it->getIntensity() = linear2log(it->getIntensity(), is_intensity_axis_percent, data_area_.maxY());
					//				cerr << it->getIntensity() << endl;
				}
			}
		}
		else
		{
			for (SpectrumIteratorType it=currentDataSet()[0].begin();it!=currentDataSet()[0].end();++it)
			{
				if (it->getIntensity()!=0)
				{
					//					cerr << it->getIntensity() << " => ";
					it->getIntensity() = log2linear(it->getIntensity(), is_intensity_axis_percent, old_max_intensity_);
					//					cerr << it->getIntensity() << endl;
				}
			}
		}
		currentDataSet()[0].updateRanges();
	}
	
	void Spectrum1DCanvas::scaleAllData_(bool is_log)
	{
		//cout << "SCALE ALL DATA"<<endl;
		//abort if there is no data
		if (getDataSetCount()==0)	return;
	
		bool is_intensity_axis_percent = !isAbsoluteIntensity();
	
		for (UnsignedInt index = 0; index < getDataSetCount(); ++index)
		{
			if (getDataSet(index)[0].size()==0) continue;
	
			if (is_log)
			{
				old_max_intensity_ = data_area_.maxY();
				for (SpectrumIteratorType it=getDataSet(index)[0].begin();it!=getDataSet(index)[0].end();++it)
				{
					if (it->getIntensity()!=0)
					{
						//cerr<<  it->getIntensity() << " => ";
						it->getIntensity() = linear2log(it->getIntensity(), is_intensity_axis_percent, data_area_.maxY());
						//cerr << it->getIntensity() << endl;
					}
				}
			}
			else
			{
				for (SpectrumIteratorType it=getDataSet(index)[0].begin();it!=getDataSet(index)[0].end();++it)
				{
					if (it->getIntensity()!=0)
					{
						//cerr << it->getIntensity() << " => ";
						it->getIntensity() = log2linear(it->getIntensity(), is_intensity_axis_percent, old_max_intensity_);
						//cerr << it->getIntensity() << endl;
					}
				}
			}
			getDataSet(index)[0].updateRanges();
			//cout << "/SCALE ALL DATA"<<endl;
		}
	
		//update visible area
		data_area_.setMinY( (is_log)? 0.0 : 1.0);
		data_area_.setMaxY(getDataSet(0)[0].getMaxInt());
	
		for (UnsignedInt index = 1; index < getDataSetCount(); ++index)
		{
			if (data_area_.maxY() < getDataSet(index)[0].getMaxInt())
			{
				data_area_.setMaxY(getDataSet(index)[0].getMaxInt());
			}
		}
	
		// extend region for cosmectical reasons -> peak with smallest/highest position won't fall on diagramm border on maximized view
		if (is_log)
		{
			//double min = log2linear(data_area_.minY(),is_intensity_axis_percent, old_max_intensity_);
			double max = log2linear(data_area_.maxY(),is_intensity_axis_percent, old_max_intensity_);
	 
			//data_area_.setMaxY( linear2log(max + 0.002*(max - min), is_intensity_axis_percent, old_max_intensity_));
			data_area_.setMaxY( linear2log(max, is_intensity_axis_percent, old_max_intensity_));
		}else
		{
			//data_area_.setMaxY( data_area_.maxY() + 0.002*(data_area_.maxY()  - data_area_.minY()));
			data_area_.setMaxY( data_area_.maxY());
		}
	
		visible_area_.setMaxY(data_area_.maxY());
		visible_area_.setMinY(data_area_.minY());
													
		//set displayed intensity range
		min_disp_ints_[current_data_] = 0;
		max_disp_ints_[current_data_] = data_area_.maxY();
	}
	
	void Spectrum1DCanvas::intensityAxisAbsolute()
	{
		if (!isAbsoluteIntensity())
		{
			// Rescale data before switching to absolute log scale
			if (intensity_modification_ == IM_LOG)
			{
				scaleData_(false);
	 			intensityModificationChange_();
			}
			absolute_intensity_ = true;
			invalidate_();
		}
	}
	
	void Spectrum1DCanvas::intensityAxisRelative()
	{
		if (isAbsoluteIntensity())
		{
			// Rescale data before switching to absolute log scale
			if (intensity_modification_ == IM_LOG)
			{
				scaleData_(false); 
				intensityModificationChange_();
			}
			absolute_intensity_ = false;
			invalidate_();
		}
	}
	
	bool Spectrum1DCanvas::isAbsoluteIntensity() const
	{
		return absolute_intensity_;
	}
	
	
	
	void Spectrum1DCanvas::changeVisibleArea_(const Spectrum1DCanvas::AreaType& new_area)
	{
		//prevent deadlock
		if (new_area==visible_area_)
		{
			return;
		}
		
		//store old zoom state
		if (zoom_timeout_)
		{
			zoom_stack_.push(visible_area_);
		}
		
		visible_area_ = new_area;
		
		updateScrollbars_();
		
		zoom_timeout_ = false;
		QTimer::singleShot(2000, this, SLOT(timeoutZoom_()));
	
		setBounds_();
		emit visibleAreaChanged(new_area);
		recalculate_ = true;
		invalidate_();
	}
	
	// destructor
	Spectrum1DCanvas::~Spectrum1DCanvas()
	{
		
	}
	
	std::vector<Spectrum1DCanvas::SpectrumIteratorType> Spectrum1DCanvas::getSelectedPeaks()
	{
		std::vector<SpectrumIteratorType> result = selected_peaks_;
		
		//to also have the last peak of the spectrum as border: add the peak BEFORE currentDataSet()[0].end()
		if (!currentDataSet()[0].empty())
		{
			result.push_back((currentDataSet()[0].end() - 1));
		}
	
		return result; 
	}
	
	bool Spectrum1DCanvas::getSnapToMax()
	{
		return snap_to_max_mode_;
	}
	
	void Spectrum1DCanvas::setSnapToMax(bool b)
	{
		snap_to_max_mode_ = b;
		setBounds_();
		emit visibleAreaChanged(getDataArea_());
		invalidate_();
	}
	
	double Spectrum1DCanvas::getSnapFactor()
	{
		return snap_factor_;
	}
	
	PreferencesDialogPage* Spectrum1DCanvas::createPreferences(QWidget* parent)
	{
		return new Spectrum1DCanvasPDP(this, parent);
	}
	
	void Spectrum1DCanvas::clearHighlighting()
	{
		if (is_highlighted_) {
			is_highlighted_ = false;
			emit sendCursorStatus();
			nearest_peak_ = currentDataSet()[0].end();
			invalidate_();
		}
	}

	SignedInt Spectrum1DCanvas::finishAdding()
	{
		UnsignedInt index = getDataSetCount()-1;
		current_data_ = index;
		currentDataSet().updateRanges();
		
		if (currentDataSet().size()==0 || currentDataSet().getSize()==0)
		{
			datasets_.resize(getDataSetCount()-1);
			current_data_ = index-1;
			return -1;
		}
		
		//add display intensity ranges
		min_disp_ints_.resize(index+1);
		max_disp_ints_.resize(index+1);
	
		//set visibility to true
		layer_visible_.push_back(true);

		
		if (getIntensityModification()==IM_LOG)
		{
			scaleData_(true);
		}
	
		// sort peaks in non descending order of position
		currentDataSet()[0].getContainer().sortByNthPosition(0);
		
		//if this is the first spectrum, set the min/max pos/int to the first value
		//update visible area
		if (index == 0)
		{
			data_area_.setMinX(currentDataSet().getMinMZ());
			data_area_.setMaxX(currentDataSet().getMaxMZ());
			data_area_.setMinY(0.0);  // minimal intensity always 0.0
			data_area_.setMaxY(currentDataSet().getMaxInt());
	
			//extend region for cosmectical reasons -> peak with smallest/highest position 
			//won't fall on diagramm border on maximized view
			data_area_.setMinX(data_area_.minX() - 0.002 * data_area_.width());
			data_area_.setMaxX(data_area_.maxX() + 0.002 * data_area_.width());
			data_area_.setMaxY(data_area_.maxY() + 0.002 * data_area_.height());
		}
	
		if (data_area_.minX() > currentDataSet().getMinMZ())
		{
			data_area_.setMinX(currentDataSet().getMinMZ());
			data_area_.setMinX(data_area_.minX() - 0.002 * data_area_.width());
		}
		if (data_area_.maxX() < currentDataSet().getMaxMZ())
		{
			data_area_.setMaxX(currentDataSet().getMaxMZ());
			data_area_.setMaxX(data_area_.maxX() + 0.002 * data_area_.width());
		}
	
		if (data_area_.maxY() < currentDataSet().getMaxInt())
		{
			data_area_.setMaxY(currentDataSet().getMaxInt());
			data_area_.setMaxY(data_area_.maxY() + 0.002 * data_area_.height());
			max_layer_ = index;
		}
	
		//set displayed intensity range
		min_disp_ints_[current_data_] = currentDataSet().getMinInt();
		max_disp_ints_[current_data_] = currentDataSet().getMaxInt();
	
		//add new values to visible_begin_ and visible_end_
		visible_begin_.push_back(currentDataSet()[0].begin());
		visible_end_.push_back(currentDataSet()[0].end());
	
		//add new draw mode
		draw_modes_.push_back(DM_PEAKS);
	
		changeVisibleArea_(data_area_);
	
		//
		if (data_area_.width() < 1.0)
		{
			changeVisibleArea_(data_area_.minX() -1.0, data_area_.maxX() + 1.0);
		}
		else
		{
			changeVisibleArea_(data_area_.minX(), data_area_.maxX());
		}
		
		showGridLines(show_grid_);
	
		return index;
	}

}//Namespace




