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

// Qt
#include <qmessagebox.h>
#include <qaction.h>
#include <qtimer.h>

// OpenMS
#include <OpenMS/VISUAL/PeakIcon.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DCanvasPDP.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	using namespace Internal;
		
	Spectrum1DCanvas::Spectrum1DCanvas(QWidget* parent, const char* name, WFlags f)
		: SpectrumCanvas(parent, name, f | WRepaintNoErase)
	{
		
	}
	
	//change the current layer
	void Spectrum1DCanvas::activateLayer(int layer_index)
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()) || layer_index==int(current_layer_))
		{
			return ;
		}
		
		current_layer_ = layer_index;
			
		// no peak is selected
		nearest_peak_ = currentPeakData_()[0].end();
		selected_peaks_.clear();
		selected_peaks_.push_back(currentPeakData_()[0].begin());
		
		emit layerActivated(this);
	}
	
	void Spectrum1DCanvas::paintEvent(QPaintEvent* e)
	{
		SpectrumCanvas::paintEvent(e);
		
		QPainter painter(this);
		rubber_band_.draw(painter);
		painter.end();
	}
	
	void Spectrum1DCanvas::setVisibleArea(DRange<2> range)
	{
		changeVisibleArea_(AreaType(range.minX(), visible_area_.minY(), range.maxX(), visible_area_.maxY()));
	}
	
	void Spectrum1DCanvas::changeVisibleArea_(double lo, double hi, bool add_to_stack)
	{
		changeVisibleArea_(AreaType(lo, visible_area_.minY(), hi, visible_area_.maxY()), add_to_stack);
	}
	
	QPoint Spectrum1DCanvas::dataToWidget_(const PeakType& peak)
	{
		return SpectrumCanvas::dataToWidget_(peak.getPosition()[0], snap_factor_*percentage_factor_*peak.getIntensity());
	}
	
	//////////////////////////////////////////////////////////////////////////////////
	// Qt events
	
	void Spectrum1DCanvas::mousePressEvent( QMouseEvent *e)
	{
		// get mouse position in widget coordinates
		QPoint p = e->pos();
		last_mouse_pos_ = p;
	
		if (e->button() == LeftButton)
		{
			PointType pos = widgetToData_(p);
			if (e->state() & QMouseEvent::ShiftButton)
			{
				if (action_mode_ == AM_ZOOM)
				{
					rubber_band_.show();
					rubber_band_.setTopLeft(p);
					rubber_band_.setBottomRight(p);
					rubber_band_.updateRegion(this);
				}
			}
			else
			{
				if (action_mode_ == AM_ZOOM)
				{
					rubber_band_.show();	
					rubber_band_.setTopLeft(SpectrumCanvas::dataToWidget_(pos.X(), overall_data_range_.maxY()));
					rubber_band_.setBottomRight(SpectrumCanvas::dataToWidget_(pos.X(), overall_data_range_.minY()));
	
					rubber_band_.updateRegion(this);
				}
			}
			if (action_mode_ == AM_TRANSLATE) setCursor(cursor_translate_in_progress_);
		}
	}
	
	
	void Spectrum1DCanvas::mouseReleaseEvent(QMouseEvent *e)
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
				else if (nearest_peak_ != currentPeakData_()[0].end())
				{
					if (!nearest_peak_->metaValueExists(4) || (UnsignedInt)(nearest_peak_->getMetaValue(4)) == PeakIcon::IT_NOICON)
					{
						nearest_peak_->setMetaValue(4, SignedInt(PeakIcon::IT_ELLIPSE));
						selected_peaks_.push_back(nearest_peak_);
	
						ostringstream msg;
						msg << "Selected peak at position " << nearest_peak_->getPosition()[0]  << " (" << (selected_peaks_.size()-1);
						msg << " peaks selected altogether.)";
						emit sendStatusMessage(msg.str(), 5000);
					}
					else
					{
						nearest_peak_->setMetaValue(4,SignedInt(PeakIcon::IT_NOICON));
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
					rubber_band_.updateRegion(this);
					QRect rect = rubber_band_.getRect();
					rubber_band_.setRect(QRect(0,0,0,0));
	
					if (rect.width() > 4 && rect.height() > 4)   // free zoom
					{
						Spectrum1DCanvas::AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
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
							changeVisibleArea_(area.minX(), area.maxX(), true);
						}
					}
					else                                        // position axis only zoom
					{
						double position = widgetToData_(rect.topLeft()).X();
						double delta = (visible_area_.maxX()-visible_area_.minX())/4.0;
						changeVisibleArea_(position - delta, position + delta, true);
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
				setCursor(cursor_translate_);
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
	
	void Spectrum1DCanvas::mouseDoubleClickEvent( QMouseEvent *e)
	{
		if (e->button() == LeftButton && action_mode_ == AM_SELECT)
		{
			SpectrumIteratorType i = findPeakAtPosition(e->pos());
			if (i != currentPeakData_()[0].end())
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
	
	void Spectrum1DCanvas::mouseMoveEvent( QMouseEvent *e)
	{
		// mouse position relative to the diagram widget
		QPoint p = e->pos();
	
		switch (action_mode_)
		{
			case AM_SELECT:
			{ 
				emit sendCursorStatus();
	 			setCursor(Qt::ArrowCursor);
				nearest_peak_ = findPeakAtPosition(p);
				invalidate_();
				break;
			}
			case AM_ZOOM:
			{
				setCursor(Qt::CrossCursor);
				PointType pos = widgetToData_(p);
	
				if (e->state() & LeftButton)
				{
					if (e->state() & QMouseEvent::ShiftButton) // free zoom
					{
						rubber_band_.updateRegion(this);
						rubber_band_.setBottomRight(p);
						rubber_band_.updateRegion(this);
					}
					else // zoom on position axis only
					{
						rubber_band_.updateRegion(this);
						rubber_band_.setBottomRight(SpectrumCanvas::dataToWidget_(pos.X(), overall_data_range_.minY()));
						rubber_band_.updateRegion(this);
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
				setCursor(cursor_translate_);
				// Translation of visible area
				if(e->state() & QMouseEvent::LeftButton)
				{
					setCursor(cursor_translate_in_progress_);
					// translation in data metric
					double shift = widgetToData_(last_mouse_pos_).X() - widgetToData_(p).X();
					double newLo = visible_area_.minX() + shift;
					double newHi = visible_area_.maxX() + shift;
					// check if we are falling out of bounds
					if (newLo < overall_data_range_.minX())
					{
						newLo = overall_data_range_.minX();
						newHi = newLo + visible_area_.width();
					}
					if (newHi > overall_data_range_.maxX())
					{
						newHi = overall_data_range_.maxX();
						newLo = newHi - visible_area_.width();
					}
					//chage data area
					changeVisibleArea_(newLo, newHi);
					last_mouse_pos_=p;
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
		//		find the 2 iterators in the layer that points to the first peak that falls into this interval
		//		and the first, that falls off the interval
		//	Depending on how many peaks (if any) lie between this two iterators further decisions have to be made
		//  Case 0: no peak
		//	Case 1: one peak => we already found the nearest peak
		// 	Case 2: several peaks => as all peaks between these two iterators will be projected on the same line on the screen further
		//					decision must me made (based on P.y) (see code)
		// now a slightly modified version is used. the only difference is, that the source interval
		// has a length of 3 pixels (to simplify grabbing of a peak).
		
		// get the interval (in diagramm metric) that will be projected on screen coordinate p.x() or p.y() (depending on orientation)
		PointType lt = widgetToData_(p - QPoint(1, 1));
		PointType rb = widgetToData_(p + QPoint(1, 1));
		double interval_start = min(lt.X(),rb.X());
		double interval_end = max(lt.X(),rb.X());
		
		// debug code:
		//	cout << "Intervall start: " << interval_start << endl;
		//	cout << "Intervall end:   " << interval_end << endl;
	
		// get iterator on first peak with higher position than interval_start
		PeakType temp;
		temp.getPosition()[0] = interval_start;
		SpectrumIteratorType left_it = lower_bound(visible_begin_[current_layer_], visible_end_[current_layer_], temp, PeakType::PositionLess());
	
		// get iterator on first peak with higher position than interval_end
		temp.getPosition()[0] = interval_end;
		SpectrumIteratorType	right_it = lower_bound(visible_begin_[current_layer_], visible_end_[current_layer_], temp, PeakType::PositionLess());
	
		//debug code:
		//	cout << "left_it and *left_it: "  << left_it->getPosition()[0] << " " << left_it->getIntensity() << endl;
		//	cout << "right_it and *right_it: "  << right_it->getPosition()[0] << " " << right_it->getIntensity() << endl;
	
		if (left_it == right_it)
		{
		//	cout << "case 0: no peak" << endl;	// debug code
			return currentPeakData_()[0].end();  // both are equal => no peak falls into this interval
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
	
		double dest_interval_start = SpectrumCanvas::dataToWidget_(0, overall_data_range_.minY()).y();
		double dest_interval_end = SpectrumCanvas::dataToWidget_(0, overall_data_range_.maxY()).y();
		//double dest_interval_start, dest_interval_end;
		// select source interval start and end depending on diagram orientation
	
		if (isMzToXAxis())
		{
			// select destination interval start and end depending on diagram orientation AND axis direction
			dest_interval_start = height();
			dest_interval_end = 0;
		}
		else
		{
			// select destination interval start and end depending on diagram orientation and axis direction
			dest_interval_start = 0;
			dest_interval_end = width();
		}
	
		int nearest_intensity = static_cast<int>(intervalTransformation(nearest_it->getIntensity(), visible_area_.minY(),
		                                                                 visible_area_.maxY(), dest_interval_start, dest_interval_end));
		int current_intensity;
	
		for (SpectrumIteratorType it = left_it; it != right_it; it++)
		{
			current_intensity = static_cast<int>(intervalTransformation(it->getIntensity(), visible_area_.minY(), visible_area_.maxY(),
			                                                             dest_interval_start, dest_interval_end));
			if (isMzToXAxis())
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
	
	void Spectrum1DCanvas::removeLayer(int layer_index)
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()))
		{
			return;
		}
	
		//remove settings
		layers_.erase(layers_.begin()+layer_index);
		draw_modes_.erase(draw_modes_.begin()+layer_index);
	
		//refresh values of visible_begin_ and visible_end_
		visible_begin_.clear();
		visible_end_.clear();
		for (UnsignedInt index=0; index < getLayerCount(); ++index)
		{
			visible_begin_.push_back(getPeakData_(index)[0].begin());
			visible_end_.push_back(getPeakData_(index)[0].end());
		}
	
		//update current layer
		if (current_layer_ >= getLayerCount())
		{
			current_layer_ = getLayerCount()-1;
		}
	
		if (layers_.empty())
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
		
		invalidate_();
	}

	void Spectrum1DCanvas::setDrawMode(DrawModes mode)
	{
		if (draw_modes_[current_layer_]!=mode)
		{
			draw_modes_[current_layer_] = mode;
			invalidate_();
		}
	}

	Spectrum1DCanvas::DrawModes Spectrum1DCanvas::getDrawMode() const
	{ 
		return draw_modes_[current_layer_]; 
	}

	
	//////////////////////////////////////////////////////////////////////////////////
	// data plotting
	
	void Spectrum1DCanvas::drawIcon(const PeakType& peak, const QPoint& p)
	{
		// draw icon associated to current peak
		if (peak.metaValueExists(4))
		{
			painter_.setPen(icon_pen_);	
			PeakIcon::drawIcon((PeakIcon::Icon)(UnsignedInt)(peak.getMetaValue(4)),painter_,QRect(p.x() - 5, p.y() - 5, 10, 10));
			painter_.setPen(norm_pen_);
		}
	}
	
	void Spectrum1DCanvas::drawPeaks_(UnsignedInt index)
	{
		painter_.save();
		painter_.setBrush(NoBrush);
		painter_.setPen(norm_pen_);
		
		//Factor to stretch the log value to the shown intensity interval
		float log_factor = getPeakData(index).getMaxInt()/log(getPeakData(index).getMaxInt());
		
		QPoint p, p0;
		bool custom_color;
		double min_int = getLayer(index).min_int;
		double max_int = getLayer(index).max_int;
		for (SpectrumIteratorType it = visible_begin_[index]; it != visible_end_[index]; ++it)
		{
			if (it->getIntensity() < min_int || it->getIntensity() > max_int)
			{
				continue;
			}
			if (intensity_mode_==IM_LOG)
			{
				p = SpectrumCanvas::dataToWidget_(it->getPosition()[0], log(it->getIntensity()+1)*log_factor);
			}
			else
			{
				p = dataToWidget_(*it);
			}
			p0 = SpectrumCanvas::dataToWidget_(it->getPosition()[0], 0.0f);
			
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
		float log_factor = getPeakData(index).getMaxInt()/log(getPeakData(index).getMaxInt());
		
		//cases where 1 or 0 points are shown
		if (visible_begin_[index]==visible_end_[index])
		{
			// check cases where no peak at all is visible
			if (visible_begin_[index] == getPeakData_(index)[0].end()) 
			{
				painter_.restore();
				return;
			}
			if (visible_end_[index] == getPeakData_(index)[0].begin())
			{
				painter_.restore();
				return;
			}
			// draw line (clipping performed by Qt on both sides)
			painter_.drawLine(dataToWidget_(*(visible_begin_[index] - 1)), dataToWidget_(*visible_begin_[index]));
			return;
		}
	
		// connect peaks in visible area; (no clipping needed)
		bool firstPoint=true;
		QPoint p;
		bool custom_color;
		for (SpectrumIteratorType it = visible_begin_[index]; it != visible_end_[index]; it++)
		{
			if (intensity_mode_==IM_LOG)
			{
				p = SpectrumCanvas::dataToWidget_(it->getPosition()[0], log(it->getIntensity()+1)*log_factor);
			}
			else
			{
				p = dataToWidget_(*it);
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
		if (visible_begin_[index] > getPeakData_(index)[0].begin())
		{
			painter_.drawLine(dataToWidget_(*(visible_begin_[index]-1)), dataToWidget_(*(visible_begin_[index])));
		}
	
		// clipping on right side
		if (visible_end_[index] < getPeakData_(index)[0].end())
		{
			painter_.drawLine( dataToWidget_(*(visible_end_[index]-1)), dataToWidget_(*(visible_end_[index])));
		}
	
		painter_.restore();
	}
	
	void Spectrum1DCanvas::setMainPreferences(const Param& prefs)
	{
		SpectrumCanvas::setMainPreferences(prefs);
		if (getPrefAsString("Preferences:1D:Mapping:MappingOfMzTo") != "X-Axis")
		{
			mzToXAxis(false);
		}
	}
	
	void Spectrum1DCanvas::invalidate_()
	{
		// get color settings
		setPaletteBackgroundColor(QColor(getPrefAsString("Preferences:1D:BackgroundColor").c_str()));
		setBackgroundColor(QColor(getPrefAsString("Preferences:1D:BackgroundColor").c_str()));
	
		norm_pen_ = QPen(QColor(getPrefAsString("Preferences:1D:PeakColor").c_str()), pen_width_);
		high_pen_ = QPen(QColor(getPrefAsString("Preferences:1D:HighColor").c_str()), pen_width_+2);
		icon_pen_ = QPen(QColor(getPrefAsString("Preferences:1D:IconColor").c_str()), pen_width_);
	
		//repaint content
		painter_.begin(buffer_);
		// flushing buffer_ with widget color
		painter_.fillRect(0, 0, width(), height(), backgroundColor());
		
		emit recalculateAxes();
		paintGridLines_(&painter_);
		
		for (UnsignedInt i=0; i< getLayerCount();++i)
		{
			if (!getLayer(i).visible)
			{
				continue;
			}
	
			if (intensity_mode_ == IM_PERCENTAGE)
			{
				percentage_factor_ = overall_data_range_.max()[1]/getPeakData(i)[0].getMaxInt();
			}
			else 
			{
				percentage_factor_ = 1.0;
			}
			
			switch (draw_modes_[i])
			{
				case DM_PEAKS:
					drawPeaks_(i);
					break;
				case DM_CONNECTEDLINES:
					drawConnectedLines_(i);
					break;
				default:
					throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}
		}
	
		painter_.end();
		repaint(false);
	}
	
	void Spectrum1DCanvas::changeVisibleArea_(const AreaType& new_area, bool add_to_stack)
	{
		//prevent deadlock
		if (new_area==visible_area_)
		{
			return;
		}
	
		if (!layers_.empty())
		{
			// get iterators on peaks that outline the visible area
			for (UnsignedInt i=0; i<getLayerCount();++i)
			{
				visible_begin_[i] = getPeakData_(i)[0].MZBegin(new_area.minX());
				visible_end_[i]   = getPeakData_(i)[0].MZBegin(new_area.maxX());
			}
			
			if (action_mode_ != AM_SELECT)
			{
				nearest_peak_ = visible_end_[current_layer_];
			}

			recalculateSnapFactor_();	
		}

		SpectrumCanvas::changeVisibleArea_(new_area, add_to_stack);
	}
	
	// destructor
	Spectrum1DCanvas::~Spectrum1DCanvas()
	{
		
	}
	
	vector<Spectrum1DCanvas::SpectrumIteratorType> Spectrum1DCanvas::getSelectedPeaks()
	{
		vector<SpectrumIteratorType> result = selected_peaks_;
		
		//to also have the last peak of the spectrum as border: add the peak BEFORE getCurrentPeakData()[0].end()
		if (!getCurrentPeakData()[0].empty())
		{
			result.push_back((currentPeakData_()[0].end() - 1));
		}
	
		return result; 
	}
	
	PreferencesDialogPage* Spectrum1DCanvas::createPreferences(QWidget* parent)
	{
		return new Spectrum1DCanvasPDP(this, parent);
	}

	SignedInt Spectrum1DCanvas::finishAdding(float low_intensity_cutoff)
	{
		current_layer_ = getLayerCount()-1;
		currentPeakData_().updateRanges();
		
		if (getCurrentPeakData().size()==0 || getCurrentPeakData().getSize()==0)
		{
			layers_.resize(getLayerCount()-1);
			current_layer_ = current_layer_-1;
			return -1;
		}

		//set displayed intensity range
		getCurrentLayer_().min_int = low_intensity_cutoff;
		getCurrentLayer_().max_int = getCurrentPeakData().getMaxInt();
	
		//add new values to visible_begin_ and visible_end_
		visible_begin_.push_back(currentPeakData_()[0].begin());
		visible_end_.push_back(currentPeakData_()[0].end());
	
		//add new draw mode
		draw_modes_.push_back(DM_PEAKS);
		//estimate peak type
		PeakTypeEstimator pte;
		if (pte.estimateType(currentPeakData_()[0].begin(),currentPeakData_()[0].end()) == SpectrumSettings::RAWDATA)
		{
			draw_modes_.back() = DM_CONNECTEDLINES;
		}
	
		// sort peaks in accending order of position
		currentPeakData_()[0].getContainer().sortByPosition();
		
		//update ranges
		recalculateRanges_(0,2,1);
		overall_data_range_.setMinY(0.0);  // minimal intensity always 0.0
		float width = overall_data_range_.width();
		overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
		overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
		overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());
		
		resetZoom();
		
		emit layerActivated(this);
		
		return current_layer_;
	}

  void Spectrum1DCanvas::recalculateSnapFactor_()
  {
		if (intensity_mode_ == IM_SNAP) 
		{
			double local_max  = -numeric_limits<double>::max();
			for (UnsignedInt i=0; i<getLayerCount();++i)
			{
				SpectrumIteratorType tmp  = max_element(visible_begin_[i], visible_end_[i], PeakType::IntensityLess());
				if (tmp->getIntensity() > local_max) 
				{
					local_max = tmp->getIntensity();
				}
			}
			snap_factor_ = overall_data_range_.max()[1]/local_max;
		}
		else
		{ 
			snap_factor_ = 1.0;
		}  	
  }

	void Spectrum1DCanvas::updateScrollbars_()
	{
		if (isMzToXAxis())
		{
			emit updateHScrollbar(overall_data_range_.min()[0],visible_area_.min()[0],visible_area_.max()[0],overall_data_range_.max()[0]);
			emit updateVScrollbar(1,1,1,1);
		}
		else
		{
			emit updateHScrollbar(1,1,1,1);
			emit updateVScrollbar(overall_data_range_.min()[0],visible_area_.min()[0],visible_area_.max()[0],overall_data_range_.max()[0]);
		}
	}

	void Spectrum1DCanvas::horizontalScrollBarChange(int value)
	{
		if (isMzToXAxis())
		{
			changeVisibleArea_(value, value + (visible_area_.max()[0] - visible_area_.min()[0]));
		}
	}

	void Spectrum1DCanvas::verticalScrollBarChange(int value)
	{
		if (!isMzToXAxis())
		{
			double range = (overall_data_range_.maxX() - overall_data_range_.minX())- (visible_area_.maxX() - visible_area_.minX());
			double newval = (1.0 - (double(value) - overall_data_range_.minX()) / range )* range + overall_data_range_.minX();
			//cout << value << " " << newval << " " << newval + (visible_area_.maxX() - visible_area_.minX()) << endl;
			//cout << "Min: " <<  overall_data_range_.minX() << " Range: " << range << endl << endl;
			changeVisibleArea_(newval, newval + (visible_area_.maxX() - visible_area_.minX()));
		}
	}

}//Namespace




