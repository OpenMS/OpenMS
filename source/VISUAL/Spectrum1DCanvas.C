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
		absolute_intensity_(false),
		layer_factor_(1.0),
		snap_to_max_mode_(false),
		snap_factor_(1.0)
	{
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
	
	void Spectrum1DCanvas::zoomIn(double position)
	{
		double delta, newLo, newHi;
	
		delta = (visible_area_.maxX()-visible_area_.minX())/4.0;
		newLo = position - delta;
		newHi = position + delta;
		
		changeVisibleArea_(newLo, newHi);
	}
	
	
	void Spectrum1DCanvas::zoomOut(double position)
	{
		double delta;
		double newLo;
		double newHi;
	
		delta = (visible_area_.maxX()-visible_area_.minX());
		newLo = position-delta;
		newHi = position+delta;
	
		// fit too data set bounds
		if (newLo < overall_data_range_.minX()) newLo = overall_data_range_.minX();
		if (newHi > overall_data_range_.maxX()) newHi = overall_data_range_.maxX();
	
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
	
		if (newLo < overall_data_range_.minX())
		{
			changeVisibleArea_(overall_data_range_.minX(), overall_data_range_.minX() + 2 * delta);
			return;
		}
	
		if (newHi >overall_data_range_.maxX())
		{
			changeVisibleArea_(overall_data_range_.maxX() - 2 * delta, overall_data_range_.maxX());
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
				}
			}
			else
			{
				if (action_mode_ == AM_ZOOM)
				{
					rubber_band_.show();	
					rubber_band_.setTopLeft(SpectrumCanvas::chartToWidget_(PointType(pos.X(), overall_data_range_.maxY())));
					rubber_band_.setBottomRight(SpectrumCanvas::chartToWidget_(PointType(pos.X(), overall_data_range_.minY())));
	
					rubber_band_.updateRegion(viewport());
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
	
						ostringstream msg;
						msg << "Selected peak at position " << nearest_peak_->getPosition()[0]  << " (" << (selected_peaks_.size()-1);
						msg << " peaks selected altogether.)";
						emit sendStatusMessage(msg.str(), 5000);
					}
					else
					{
						nearest_peak_->setMetaValue(4,SignedInt(Spectrum1DCanvas::IT_NOICON));
						vector<SpectrumIteratorType>::iterator it_tmp = std::find(selected_peaks_.begin(), selected_peaks_.end(), nearest_peak_);
	
						if(it_tmp != selected_peaks_.end())
						{
							selected_peaks_.erase(it_tmp);
	
							ostringstream msg;
							msg << "Deselected peak at position " << nearest_peak_->getPosition()[0]  << " (" << (selected_peaks_.size()-1);
							msg << " peaks selected altogether.)";
							emit sendStatusMessage(msg.str(), 5000);
						}
	
						//cout << "selected_peaks_.size(): " << selected_peaks_.size() << endl;
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
					rubber_band_.hide();
					rubber_band_.updateRegion(viewport());
					QRect rect = rubber_band_.getRect();
					rubber_band_.setRect(QRect(0,0,0,0));
	
					if (rect.width() > 4 && rect.height() > 4)   // free zoom
					{
						Spectrum1DCanvas::AreaType area(widgetToChart_(rect.topLeft()), widgetToChart_(rect.bottomRight()));
						if (e->state() & QMouseEvent::ShiftButton)
						{
							//check if selected area is outside visible area and correct errors
							if (area.minX() < visible_area_.minX())
							{
								area.setMinX(visible_area_.minX());
							}
							if (area.minY() < visible_area_.minY())
							{
								area.setMinY(visible_area_.minY());
							}
							if (area.maxX() > visible_area_.maxX())
							{
								area.setMaxX(visible_area_.maxX());
							}
							if (area.maxY() > visible_area_.maxY())
							{
								area.setMaxY(visible_area_.maxY());
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
						zoomIn(widgetToChart_(rect.topLeft()).X());
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
			case AM_MEASURE:
			{
				
			}	
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
				nearest_peak_ = findPeakAtPosition(p);
				invalidate_();
				break;
			}
			case AM_ZOOM:
			{
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
						rubber_band_.setBottomRight(SpectrumCanvas::chartToWidget_(PointType(pos.X(), overall_data_range_.minY())));
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
						if (newLo < overall_data_range_.minX())
						{
							changeVisibleArea_(overall_data_range_.minX(), overall_data_range_.minX()+(visible_area_.maxX()-visible_area_.minX()));
						}
						else if (newHi > overall_data_range_.maxX())
						{
							changeVisibleArea_(overall_data_range_.maxX()-(visible_area_.maxX()-visible_area_.minX()), overall_data_range_.maxX());
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
			case AM_MEASURE:
			{
				
			}
		}
	}
	
	Spectrum1DCanvas::SpectrumIteratorType Spectrum1DCanvas::findPeakAtPosition(QPoint p)
	{
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
		
		// get the interval (in diagramm metric) that will be projected on screen coordinate p.x() or p.y() (depending on orientation)
		PointType lt = widgetToChart_(p - QPoint(1, 1));
		PointType rb = widgetToChart_(p + QPoint(1, 1));
		double interval_start = min(lt.X(),rb.X());
		double interval_end = max(lt.X(),rb.X());
		
		// debug code:
		//	cout << "Intervall start: " << interval_start << endl;
		//	cout << "Intervall end:   " << interval_end << endl;
	
		// get iterator on first peak with higher position than interval_start
		PeakType temp;
		temp.getPosition()[0] = interval_start;
		SpectrumIteratorType left_it = lower_bound(visible_begin_[current_data_], visible_end_[current_data_], temp, PeakType::PositionLess());
	
		// get iterator on first peak with higher position than interval_end
		temp.getPosition()[0] = interval_end;
		SpectrumIteratorType	right_it = lower_bound(visible_begin_[current_data_], visible_end_[current_data_], temp, PeakType::PositionLess());
	
		//debug code:
		//	cout << "left_it and *left_it: "  << left_it->getPosition()[0] << " " << left_it->getIntensity() << endl;
		//	cout << "right_it and *right_it: "  << right_it->getPosition()[0] << " " << right_it->getIntensity() << endl;
	
		if (left_it == right_it)
		{
		//	cout << "case 0: no peak" << endl;	// debug code
			return currentDataSet()[0].end();  // both are equal => no peak falls into this interval
		}
	
		if (left_it == right_it-1 )
		{
		//	cout << "case 1: one peak falls into the interval " << endl; // debug code
			return left_it;
		}
	
		// debug code
		//	cout << " case 3: several peaks in the same interval" << endl;
		//	cout << left_it->getIntensity() << " " << right_it->getIntensity() << endl;
	
		SpectrumIteratorType nearest_it = left_it;
	
		double dest_interval_start = SpectrumCanvas::chartToWidget_(DPosition<2>(0, overall_data_range_.minY())).y();
		double dest_interval_end = SpectrumCanvas::chartToWidget_(DPosition<2>(0, overall_data_range_.maxY())).y();
		//double dest_interval_start, dest_interval_end;
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
	
		//remove settings
		datasets_.erase(datasets_.begin()+data_set);
		disp_ints_.erase(disp_ints_.begin()+data_set);	
		layer_visible_.erase(layer_visible_.begin()+data_set);
		draw_modes_.erase(draw_modes_.begin()+data_set);
	
		//refresh values of visible_begin_ and visible_end_
		visible_begin_.clear();
		visible_end_.clear();
		for (UnsignedInt index=0; index < getDataSetCount(); ++index)
		{
			visible_begin_.push_back(getDataSet(index)[0].begin());
			visible_end_.push_back(getDataSet(index)[0].end());
		}
	
		//update current data set
		if (current_data_ >= getDataSetCount())
		{
			current_data_ = getDataSetCount()-1;
		}
	
		if (datasets_.empty())
		{
			resetRanges_();
			return;
		}
	
		//update range area
		recalculateRanges_(0,2,1);
		overall_data_range_.setMinY(0.0);  // minimal intensity always 0.0
		float width = overall_data_range_.width();
		overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
		overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
		overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());
		
		//cout << overall_data_range_ << endl;
		
		AreaType tmp;
		tmp.assign(overall_data_range_);
		changeVisibleArea_(tmp);
	
		//
		if (overall_data_range_.maxX() - overall_data_range_.minX() <1.0)
		{
			changeVisibleArea_(overall_data_range_.minX() -1.0, overall_data_range_.maxX() + 1.0);
		}
		else
		{
			changeVisibleArea_(overall_data_range_.minX(), overall_data_range_.maxX());
		}
		showGridLines(show_grid_);
	}
	
	void Spectrum1DCanvas::setDrawMode(QAction* a)
	{
		QString name = a->name();
	
		if (name == "setPeakMode")
		{
			draw_modes_[current_data_] = DM_PEAKS;
		}
		else if (name == "setConnectedLineMode")
		{
			draw_modes_[current_data_]= DM_CONNECTEDLINES;
		}
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
		painter_.save();
		painter_.setBrush(NoBrush);
		painter_.setPen(norm_pen_);
		
		//Factor to stretch the log value to the shown intensity interval
		float log_factor = getDataSet(index).getMaxInt()/log(getDataSet(index).getMaxInt());
		
		QPoint p, p0;
		bool custom_color;
		for (SpectrumIteratorType it = visible_begin_[index]; it != visible_end_[index]; ++it)
		{
			if (it->getIntensity() < disp_ints_[index].first || it->getIntensity() > disp_ints_[index].second)
			{
				continue;
			}
			if (intensity_modification_==IM_NONE)
			{
				p = chartToWidget_(*it);
			}
			else
			{
				p = SpectrumCanvas::chartToWidget_(PointType(it->getPosition()[0], log(it->getIntensity()+1)*log_factor));
			}
			p0 = SpectrumCanvas::chartToWidget_(PointType(it->getPosition()[0], 0.0f));
			
			// highlight selected peak
			if (it->getPosition()[0] == nearest_peak_->getPosition()[0])
			{
				painter_.setPen(high_pen_);
				painter_.drawLine(p0, p);
				painter_.setPen(norm_pen_);
				emit sendCursorStatus( it->getPosition()[0], it->getIntensity());
			}
			else
			{
				// custom peak color
				custom_color = it->metaValueExists(5);
				if (custom_color)
				{
					QPen pen(QColor(string(it->getMetaValue(5)).c_str()), pen_width_);
					painter_.setPen(pen);
				}
				painter_.drawLine(p0, p);
				if (custom_color)
				{
					painter_.setPen(norm_pen_);
				}
			}
			//draw icon if necessary
			drawIcon(*it, p);
		}
		painter_.restore();
	}
	
	void Spectrum1DCanvas::drawConnectedLines_(UnsignedInt index)
	{
		painter_.save();
		painter_.setBrush(NoBrush);
		painter_.setPen(norm_pen_);

		//Factor to stretch the log value to the shown intensity interval
		float log_factor = getDataSet(index).getMaxInt()/log(getDataSet(index).getMaxInt());
		
		//cases where 1 or 0 points are shown
		if (visible_begin_[index]==visible_end_[index])
		{
			// check cases where no peak at all is visible
			if (visible_begin_[index] == getDataSet(index)[0].end()) 
			{
				painter_.restore();
				return;
			}
			if (visible_end_[index] == getDataSet(index)[0].begin())
			{
				painter_.restore();
				return;
			}
			// draw line (clipping performed by Qt on both sides)
			painter_.drawLine(chartToWidget_(*(visible_begin_[index] - 1)), chartToWidget_(*visible_begin_[index]));
			return;
		}
	
		// connect peaks in visible area; (no clipping needed)
		bool firstPoint=true;
		QPoint p;
		bool custom_color;
		for (SpectrumIteratorType it = visible_begin_[index]; it != visible_end_[index]; it++)
		{
			if (intensity_modification_==IM_NONE)
			{
				p = chartToWidget_(*it);
			}
			else
			{
				p = SpectrumCanvas::chartToWidget_(PointType(it->getPosition()[0], log(it->getIntensity()+1)*log_factor));
			}

			// connect lines
			if (firstPoint)
			{
				
				painter_.moveTo(p);
				firstPoint = false;
			} 
			else
			{
				// custom peak color
				custom_color = it->metaValueExists(5);
				if (custom_color)
				{
					QPen pen(QColor(string(it->getMetaValue(5)).c_str()), pen_width_);
					painter_.setPen(pen);
				}
				painter_.lineTo(p);
				if (custom_color)
				{
					painter_.setPen(norm_pen_);
				}
			};
	
			// highlight selected peak
			if (it->getPosition()[0] == nearest_peak_->getPosition()[0])
			{
				painter_.save();
				painter_.setPen(high_pen_);
				painter_.drawLine(p.x(), p.y()-4, p.x(), p.y()+4);
				painter_.drawLine(p.x()-4, p.y(), p.x()+4, p.y());
				painter_.restore();

				emit sendCursorStatus( it->getPosition()[0], it->getIntensity());
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
	
	void Spectrum1DCanvas::updateVisibleAreaBounds_()
	{
		if (!datasets_.empty())
		{
			// get iterators on peaks that outline the visible area
			for (UnsignedInt i=0; i<getDataSetCount();++i)
			{
				visible_begin_[i] = getDataSet(i)[0].MZBegin(visible_area_.minX());
				visible_end_[i]   = getDataSet(i)[0].MZBegin(visible_area_.maxX());
			}
	
			// If snap-to-max-mode is on: find local max and set factor to increase all data appropriately
			if (snap_to_max_mode_) 
			{
				double local_max  = -numeric_limits<double>::max();
				for (UnsignedInt i=0; i<getDataSetCount();++i)
				{
					SpectrumIteratorType tmp  = max_element(visible_begin_[i], visible_end_[i], PeakType::IntensityLess());
					if (tmp->getIntensity() > local_max) 
					{
						local_max = tmp->getIntensity();
					}
				}
				snap_factor_ = 1.0*overall_data_range_.max()[1]/local_max;
			}
			else
			{ 
				snap_factor_ = 1.0;
			}
			
			if (action_mode_ != AM_SELECT)
			{
				nearest_peak_ = visible_end_[current_data_];
			}
		}
	}
	
	void Spectrum1DCanvas::invalidate_()
	{
		//cout << "INVALIDATE"<<endl;
		updateVisibleAreaBounds_();
	
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
				layer_factor_ = 1.0*overall_data_range_.max()[1]/getDataSet(i)[0].getMaxInt();
	
			switch (draw_modes_[i])
			{
				case DM_PEAKS:
					drawPeaks_(i);
				break;
				case DM_CONNECTEDLINES:
					drawConnectedLines_(i);
				break;
			}
		}
	
		painter_.end();
		viewport()->repaint(false);
		//cout << "/INVALIDATE"<<endl;
	}
	
	
	void Spectrum1DCanvas::intensityModificationChange_()
	{
		invalidate_();
	}
	
	void Spectrum1DCanvas::intensityAxisAbsolute()
	{
		if (!isAbsoluteIntensity())
		{
			absolute_intensity_ = true;
			invalidate_();
		}
	}
	
	void Spectrum1DCanvas::intensityAxisRelative()
	{
		if (isAbsoluteIntensity())
		{
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
	
		updateVisibleAreaBounds_();
		emit visibleAreaChanged(new_area);
		recalculate_ = true;
		invalidate_();
	}
	
	// destructor
	Spectrum1DCanvas::~Spectrum1DCanvas()
	{
		
	}
	
	vector<Spectrum1DCanvas::SpectrumIteratorType> Spectrum1DCanvas::getSelectedPeaks()
	{
		vector<SpectrumIteratorType> result = selected_peaks_;
		
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
		updateVisibleAreaBounds_();
		AreaType tmp;
		tmp.assign(overall_data_range_);
		emit visibleAreaChanged(tmp);
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

	SignedInt Spectrum1DCanvas::finishAdding()
	{
		current_data_ = getDataSetCount()-1;
		currentDataSet().updateRanges();
		
		if (currentDataSet().size()==0 || currentDataSet().getSize()==0)
		{
			datasets_.resize(getDataSetCount()-1);
			current_data_ = current_data_-1;
			return -1;
		}

		//set displayed intensity range
		disp_ints_.push_back(pair<float,float>(currentDataSet().getMinInt(),currentDataSet().getMaxInt()));
	
		//add new values to visible_begin_ and visible_end_
		visible_begin_.push_back(currentDataSet()[0].begin());
		visible_end_.push_back(currentDataSet()[0].end());
	
		//add new draw mode
		draw_modes_.push_back(DM_PEAKS);
	
		//set visibility to true
		layer_visible_.push_back(true);
	
		// sort peaks in accending order of position
		currentDataSet()[0].getContainer().sortByNthPosition(0);
		
		//update ranges
		recalculateRanges_(0,2,1);
		overall_data_range_.setMinY(0.0);  // minimal intensity always 0.0
		float width = overall_data_range_.width();
		overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
		overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
		overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());
		
		//cout << overall_data_range_ << endl;
		AreaType tmp;
		tmp.assign(overall_data_range_);
		changeVisibleArea_(tmp);
	
		//
		if (overall_data_range_.width() < 1.0)
		{
			changeVisibleArea_(overall_data_range_.minX() -1.0, overall_data_range_.maxX() + 1.0);
		}
		else
		{
			changeVisibleArea_(overall_data_range_.minX(), overall_data_range_.maxX());
		}
		
		showGridLines(show_grid_);
	
		return current_data_;
	}

}//Namespace




