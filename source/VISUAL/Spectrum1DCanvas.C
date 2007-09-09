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

// Qt
#include <QtGui/QMouseEvent>
#include <QtGui/QMessageBox>
#include <QtGui/QPainterPath>
#include <QtGui/QPainter>
#include <QtCore/QTime>
 
// OpenMS
#include <OpenMS/VISUAL/PeakIcon.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DPrefDialog.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/CONCEPT/TimeStamp.h>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	using namespace Internal;
		
	Spectrum1DCanvas::Spectrum1DCanvas(const Param& preferences, QWidget* parent)
		: SpectrumCanvas(preferences, parent)
	{

    //Paramater handling
    defaults_.setValue("highlighted_peak_color", "#ff0000", "Highlighted peak color.");
    defaults_.setValue("icon_color", "#000000", "Peak icon color.");
    defaults_.setValue("peak_color", "#0000ff", "Peak color.");
    defaults_.setValue("background_color", "#ffffff", "Background color.");
		defaultsToParam_();
		setName("Spectrum1DCanvas");
		setParameters(preferences);
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
		selected_peak_ = currentPeakData_()[0].end();
		selected_peaks_.clear();
		selected_peaks_.push_back(currentPeakData_()[0].begin());
		
		emit layerActivated(this);
	}
	
	void Spectrum1DCanvas::setVisibleArea(DRange<2> range)
	{
		changeVisibleArea_(AreaType(range.minX(), visible_area_.minY(), range.maxX(), visible_area_.maxY()));
	}
	
	void Spectrum1DCanvas::changeVisibleArea_(double lo, double hi, bool add_to_stack)
	{
		changeVisibleArea_(AreaType(lo, visible_area_.minY(), hi, visible_area_.maxY()), add_to_stack);
	}
	
	void Spectrum1DCanvas::dataToWidget_(const PeakType& peak, QPoint& point)
	{
		SpectrumCanvas::dataToWidget_(peak.getMZ(), snap_factor_*percentage_factor_*peak.getIntensity(), point);
	}
	
	//////////////////////////////////////////////////////////////////////////////////
	// Qt events
	
	void Spectrum1DCanvas::mousePressEvent( QMouseEvent* e)
	{
		// get mouse position in widget coordinates
		last_mouse_pos_ = e->pos();
	
		if (e->button() == Qt::LeftButton)
		{
			if (e->modifiers() & Qt::ControlModifier) //translate
			{
      	setCursor(Qt::ClosedHandCursor);
			}
			else //zoom
			{
				rubber_band_.setGeometry(e->pos().x(),e->pos().y(),0,0);
				rubber_band_.show();
			}
		}
	}

	void Spectrum1DCanvas::mouseMoveEvent(QMouseEvent* e)
	{
		// mouse position relative to the diagram widget
		QPoint p = e->pos();
	
		switch (action_mode_)
		{
			case AM_SELECT:
			{ 
				emit sendCursorStatus();
				selected_peak_ = findPeakAtPosition_(p);
				update_(__PRETTY_FUNCTION__);
				break;
			}
			case AM_ZOOM:
			{
				if(e->buttons() & Qt::LeftButton)
				{
					if (e->modifiers() & Qt::ControlModifier) //translate
					{
						// translation in data metric
						double shift = widgetToData_(last_mouse_pos_).getX() - widgetToData_(p).getX();
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
					}
					else //zoom
					{
						PointType pos = widgetToData_(p);
			
						if (e->buttons() & Qt::LeftButton)
						{
							rubber_band_.setGeometry(last_mouse_pos_.x(), 0, p.x() - last_mouse_pos_.x(), height());
							update_(__PRETTY_FUNCTION__);
						}
						emit sendCursorStatus( pos.getX() );
					}
				}
				else
				{
					if (e->modifiers() & Qt::ControlModifier) //translate
					{
						setCursor(Qt::OpenHandCursor);
					}
					else //zoom
					{
						setCursor(Qt::CrossCursor);
					}
				}
				break;
			}
			default:
				break;
		}
	}

	
	void Spectrum1DCanvas::mouseReleaseEvent(QMouseEvent* e)
	{
		switch (action_mode_)
		{
			case AM_SELECT:
			{
				// Peak selection
				if (e->button() == Qt::LeftButton && selected_peak_ != currentPeakData_()[0].end())
				{
					if (!selected_peak_->metaValueExists(4) || (UInt)(selected_peak_->getMetaValue(4)) == PeakIcon::IT_NOICON)
					{
						selected_peak_->setMetaValue(4, Int(PeakIcon::IT_ELLIPSE));
						selected_peaks_.push_back(selected_peak_);
	
						ostringstream msg;
						msg << "Selected peak at position " << selected_peak_->getMZ()  << " (" << (selected_peaks_.size()-1);
						msg << " peaks selected altogether.)";
						emit sendStatusMessage(msg.str(), 5000);
					}
					else
					{
						selected_peak_->setMetaValue(4,Int(PeakIcon::IT_NOICON));
						vector<SpectrumIteratorType>::iterator it_tmp = std::find(selected_peaks_.begin(), selected_peaks_.end(), selected_peak_);
	
						if(it_tmp != selected_peaks_.end())
						{
							selected_peaks_.erase(it_tmp);
	
							ostringstream msg;
							msg << "Deselected peak at position " << selected_peak_->getMZ()  << " (" << (selected_peaks_.size()-1);
							msg << " peaks selected altogether.)";
							emit sendStatusMessage(msg.str(), 5000);
						}
	
						//cout << "selected_peaks_.size(): " << selected_peaks_.size() << endl;
					}
					update_buffer_ = true;
					update_(__PRETTY_FUNCTION__);
				}
				break;
			}
			case AM_ZOOM:
			{
				rubber_band_.hide();
				// zoom-in-at-position or zoom-in-to-area
				if (e->button() == Qt::LeftButton)
				{
					if (e->modifiers() & Qt::ControlModifier) //translate
					{
						setCursor(Qt::OpenHandCursor);
					}
					else //zoom
					{
						QRect rect = rubber_band_.geometry();
					
						//cout << "Canvas area (x,y)-(x1,y1): " << rect.x() << "/" << rect.y() << " - " << rect.x() + rect.width() << "/" << rect.y() + rect.height() << endl;
					
						if (rect.width()!=0 && rect.height()!=0) // probably double-click -> mouseDoubleClickEvent
						{
							AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
							changeVisibleArea_(area.minX(), area.maxX(), true);
						}
					}
				}
				break;
			}
			default:
				break;
		}
	}
	
	void Spectrum1DCanvas::mouseDoubleClickEvent(QMouseEvent* e)
	{
		if (e->button() == Qt::LeftButton && action_mode_ == AM_SELECT)
		{
			SpectrumIteratorType i = findPeakAtPosition_(e->pos());
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
	
		// left-doubleclick shows the whole spectrum
		if (e->button() == Qt::LeftButton && action_mode_ == AM_ZOOM)
		{
			resetZoom();
		}
	}

	void Spectrum1DCanvas::wheelEvent(QWheelEvent* e)
	{
		if (e->delta() > 0) // forward rotation -> zoom in
		{
			DoubleReal position = visible_area_.center().getX();
			DoubleReal half_width = (visible_area_.maxX()-visible_area_.minX())/2.0*0.9;
			changeVisibleArea_(position - half_width, position + half_width, true);
		}
		else // backward rotation -> zoom out
		{
			zoomBack_();
		}
		e->accept();
	}

	Spectrum1DCanvas::SpectrumIteratorType Spectrum1DCanvas::findPeakAtPosition_(QPoint p)
	{
		//reference to the current data
		SpectrumType& spectrum = currentPeakData_()[0];
		
		// get the interval (in diagramm metric) that will be projected on screen coordinate p.x() or p.y() (depending on orientation)
		PointType lt = widgetToData_(p - QPoint(1, 1));
		PointType rb = widgetToData_(p + QPoint(1, 1));
	
		// get iterator on first peak with higher position than interval_start
		PeakType temp;
		temp.setMZ(min(lt.getX(),rb.getX()));
		SpectrumIteratorType left_it = lower_bound(spectrum.begin(), spectrum.end(), temp, PeakType::PositionLess());
	
		// get iterator on first peak with higher position than interval_end
		temp.setMZ(max(lt.getX(),rb.getX()));
		SpectrumIteratorType	right_it = lower_bound(left_it, spectrum.end(), temp, PeakType::PositionLess());
	
	
		if (left_it == right_it) // both are equal => no peak falls into this interval
		{
			return spectrum.end();
		}
	
		if (left_it == right_it-1 )
		{
			return left_it;
		}
	
		SpectrumIteratorType nearest_it = left_it;
		
		QPoint tmp;
		SpectrumCanvas::dataToWidget_(0, overall_data_range_.minY(),tmp);
		double dest_interval_start = tmp.y();
		SpectrumCanvas::dataToWidget_(0, overall_data_range_.maxY(),tmp);
		double dest_interval_end = tmp.y();
		//double dest_interval_start, dest_interval_end;
		// select source interval start and end depending on diagram orientation
	
		// select destination interval start and end
		dest_interval_start = height();
		dest_interval_end = 0;
	
		int nearest_intensity = static_cast<int>(intervalTransformation(nearest_it->getIntensity(), visible_area_.minY(),
		                                                                 visible_area_.maxY(), dest_interval_start, dest_interval_end));
		int current_intensity;
	
		for (SpectrumIteratorType it = left_it; it != right_it; it++)
		{
			current_intensity = static_cast<int>(intervalTransformation(it->getIntensity(), visible_area_.minY(), visible_area_.maxY(),
			                                                             dest_interval_start, dest_interval_end));
			if ( abs(current_intensity - p.y()) < abs(nearest_intensity - p.y()))
			{
				nearest_intensity = current_intensity;
				nearest_it = it;
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
	
		//clear other relevant variables
		selected_peaks_.clear();

		//update current layer
		if (current_layer_!=0 && current_layer_ >= getLayerCount())
		{
			current_layer_ = getLayerCount()-1;
		}
		
		//abort if there are no layers anymore
		if (layers_.empty())
		{
			overall_data_range_ = DRange<3>::empty;
			return;
		}
		
		//update nearest peak
		selected_peak_ = currentPeakData_()[0].end();
		//update selected peaks
		selected_peaks_.push_back(currentPeakData_()[0].begin());
	
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
	
		if (overall_data_range_.maxX() - overall_data_range_.minX() <1.0)
		{
			changeVisibleArea_(overall_data_range_.minX() -1.0, overall_data_range_.maxX() + 1.0);
		}
		else
		{
			changeVisibleArea_(overall_data_range_.minX(), overall_data_range_.maxX());
		}
		
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}

	void Spectrum1DCanvas::setDrawMode(DrawModes mode)
	{
		if (draw_modes_[current_layer_]!=mode)
		{
			draw_modes_[current_layer_] = mode;
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
	}

	Spectrum1DCanvas::DrawModes Spectrum1DCanvas::getDrawMode() const
	{ 
		return draw_modes_[current_layer_]; 
	}
	
	void Spectrum1DCanvas::paintEvent(QPaintEvent* e)
	{		
#ifdef DEBUG_TOPPVIEW
		cout << "BEGIN " << __PRETTY_FUNCTION__ << endl;
	  cout << "  Visible area -- m/z: " << visible_area_.minX() << " - " << visible_area_.maxX() << " int: " << visible_area_.minY() << " - " << visible_area_.maxY() << endl;
	  cout << "  Overall area -- m/z: " << overall_data_range_.min()[0] << " - " << overall_data_range_.max()[0] << " int: " << overall_data_range_.min()[1] << " - " << overall_data_range_.max()[1] << endl; 
#endif
#ifdef TIMING_TOPPVIEW
		QTime timer;
 		timer.start();
#endif
		
		QPainter painter;
		QPoint begin, end;
		float log_factor;

		if (update_buffer_)
		{
			update_buffer_ = false;
			
			painter.begin(&buffer_);

			buffer_.fill(QColor(param_.getValue("background_color").toQString()).rgb());

			emit recalculateAxes();
			paintGridLines_(painter);
			
			SpectrumIteratorType vbegin, vend;
								
			for (UInt i=0; i< getLayerCount();++i)
			{
				if (getLayer(i).visible)
				{
					QPen icon_pen = QPen(QColor(getLayer(i).param.getValue("icon_color").toQString()), 1);
					painter.setPen(QPen(QColor(getLayer(i).param.getValue("peak_color").toQString()), 1));
					if (intensity_mode_ == IM_PERCENTAGE)
					{
						percentage_factor_ = overall_data_range_.max()[1]/getPeakData(i)[0].getMaxInt();
					}
					else 
					{
						percentage_factor_ = 1.0;
					}
					
					vbegin = getPeakData_(i)[0].MZBegin(visible_area_.minX());
					vend = getPeakData_(i)[0].MZEnd(visible_area_.maxX());

					//Factor to stretch the log value to the shown intensity interval
					log_factor = getPeakData(i).getMaxInt()/log(getPeakData(i).getMaxInt());
					
					switch (draw_modes_[i])
					{
						case DM_PEAKS:
							//-----------------------------------------DRAWING PEAKS-------------------------------------------
							
							bool custom_color;
	
							for (SpectrumIteratorType it = vbegin; it != vend; ++it)
							{
								if (it->getIntensity() >= getLayer(i).min_int && it->getIntensity() <= getLayer(i).max_int)
								{
									if (intensity_mode_==IM_LOG)
									{
										SpectrumCanvas::dataToWidget_(it->getMZ(), log(it->getIntensity()+1)*log_factor,end);
									}
									else
									{
										dataToWidget_(*it,end);
									}
									SpectrumCanvas::dataToWidget_(it->getMZ(), 0.0f, begin);
									
									// draw peak
									custom_color = it->metaValueExists(5);
									if (custom_color)
									{
										painter.save();
										painter.setPen(QColor(string(it->getMetaValue(5)).c_str()));
									}
									painter.drawLine(begin, end);
									if (custom_color)
									{
										painter.restore();
									}
									
									//draw icon if necessary
									if (it->metaValueExists(4))
									{
										painter.save();
										painter.setPen(icon_pen);	
										PeakIcon::drawIcon((PeakIcon::Icon)(UInt)(it->getMetaValue(4)),painter,QRect(end.x() - 5, end.y() - 5, 10, 10));
										painter.restore();
									}
								}
							}
							//-----------------------------------------DRAWING PEAKS END-------------------------------------------
							break;
						case DM_CONNECTEDLINES:
							{
								//-------------------------------------DRAWING CONNECTED LINES-----------------------------------------
								QPainterPath path;
							
								// connect peaks in visible area; (no clipping needed)
								bool first_point=true;
								for (SpectrumIteratorType it = vbegin; it != vend; it++)
								{
									if (intensity_mode_==IM_LOG)
									{
										SpectrumCanvas::dataToWidget_(it->getMZ(), log(it->getIntensity()+1)*log_factor,begin);
									}
									else
									{
										dataToWidget_(*it, begin);
									}
						
									// connect lines
									if (first_point)
									{
										path.moveTo(begin);
										first_point = false;
									} 
									else
									{
										path.lineTo(begin);
									}
									
									// draw associated icon
									if (it->metaValueExists(4))
									{
										painter.save();
										painter.setPen(icon_pen);											
										PeakIcon::drawIcon((PeakIcon::Icon)(UInt)(it->getMetaValue(4)),painter,QRect(begin.x() - 5, begin.y() - 5, 10, 10));
										painter.restore();
									}
								}
								painter.drawPath(path);
									
								// clipping on left side
								if (vbegin != getPeakData_(i)[0].begin())
								{
									dataToWidget_(*(vbegin-1), begin);
									dataToWidget_(*(vbegin), end);
									painter.drawLine(begin, end);
								}
							
								// clipping on right side
								if (vend != getPeakData_(i)[0].end())
								{
									dataToWidget_(*(vend-1), begin);
									dataToWidget_(*(vend), end);
									painter.drawLine(begin,end);
								}
								//-------------------------------------DRAWING CONNECTED LINES END-----------------------------------------
							}
							break;
						default:
							throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
					}
				}
			}
			painter.end();
		} //if (update_buffer)
		painter.begin(this);
		
		//draw peak data
		QVector<QRect> rects = e->region().rects();
		for (int i = 0; i < (int)rects.size(); ++i)
		{
			painter.drawPixmap(rects[i].topLeft(), buffer_, rects[i]);
		}
		//draw selected peak
		if (selected_peak_!=currentPeakData_()[0].end())
		{
			painter.setPen(QPen(QColor(param_.getValue("highlighted_peak_color").toQString()), 2));		
			if (getDrawMode() == DM_PEAKS)
			{
				if (intensity_mode_==IM_LOG)
				{
					log_factor = getCurrentPeakData().getMaxInt()/log(getCurrentPeakData().getMaxInt());
					SpectrumCanvas::dataToWidget_(selected_peak_->getMZ(), log(selected_peak_->getIntensity()+1)*log_factor,end);
				}
				if (intensity_mode_==IM_PERCENTAGE)
				{
					percentage_factor_ = overall_data_range_.max()[1]/getCurrentPeakData()[0].getMaxInt();
					dataToWidget_(*selected_peak_, end);
				}
				else
				{
					dataToWidget_(*selected_peak_,end);
				}
				SpectrumCanvas::dataToWidget_(selected_peak_->getMZ(), 0.0f, begin);
				painter.drawLine(begin, end);
			}
			else if (getDrawMode() == DM_CONNECTEDLINES)
			{
				if (intensity_mode_==IM_LOG)
				{
					log_factor = getCurrentPeakData().getMaxInt()/log(getCurrentPeakData().getMaxInt());
					SpectrumCanvas::dataToWidget_(selected_peak_->getMZ(), log(selected_peak_->getIntensity()+1)*log_factor,begin);
				}
				if (intensity_mode_==IM_PERCENTAGE)
				{
					percentage_factor_ = overall_data_range_.max()[1]/getCurrentPeakData()[0].getMaxInt();
					dataToWidget_(*selected_peak_, begin);
				}
				else
				{
					dataToWidget_(*selected_peak_, begin);
				}
				painter.drawLine(begin.x(), begin.y()-4, begin.x(), begin.y()+4);
				painter.drawLine(begin.x()-4, begin.y(), begin.x()+4, begin.y());
			}
			emit sendCursorStatus( selected_peak_->getMZ(), selected_peak_->getIntensity());
		}
		painter.end();
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
#ifdef TIMING_TOPPVIEW	
		cout << "1D PaintEvent took " << timer.elapsed() << " ms" << endl;
#endif	
	}
	
	void Spectrum1DCanvas::changeVisibleArea_(const AreaType& new_area, bool add_to_stack)
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
		recalculateSnapFactor_();
		
		emit visibleAreaChanged(new_area);
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
	}
	
	// Destructor
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

	Int Spectrum1DCanvas::finishAdding(float low_intensity_cutoff)
	{
#ifdef DEBUG_TOPPVIEW
		cout << "BEGIN " << __PRETTY_FUNCTION__ << endl;
#endif
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
		
		//update nearest peak
		selected_peak_ = currentPeakData_()[0].end();
		//update selected peaks
		selected_peaks_.push_back(currentPeakData_()[0].begin());
		
		//update ranges
		recalculateRanges_(0,2,1);
		overall_data_range_.setMinY(0.0);  // minimal intensity always 0.0
		float width = overall_data_range_.width();
		overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
		overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
		overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());
		
		resetZoom();
		
		emit layerActivated(this);

#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
		
		return current_layer_;
	}

  void Spectrum1DCanvas::recalculateSnapFactor_()
  {
  	if (intensity_mode_ == IM_SNAP) 
		{
			double local_max  = -numeric_limits<double>::max();
			for (UInt i=0; i<getLayerCount();++i)
			{
				SpectrumIteratorType tmp  = max_element(getPeakData_(i)[0].MZBegin(visible_area_.minX()), getPeakData_(i)[0].MZEnd(visible_area_.maxX()), PeakType::IntensityLess());
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
		emit updateHScrollbar(overall_data_range_.min()[0],visible_area_.min()[0],visible_area_.max()[0],overall_data_range_.max()[0]);
		emit updateVScrollbar(1,1,1,1);
	}

	void Spectrum1DCanvas::horizontalScrollBarChange(int value)
	{
		changeVisibleArea_(value, value + (visible_area_.max()[0] - visible_area_.min()[0]));
	}
	
	void Spectrum1DCanvas::showCurrentLayerPreferences()
	{
		Internal::Spectrum1DPrefDialog dlg(this);
		
		ColorSelector* peak_color = dlg.findChild<ColorSelector*>("peak_color");
		ColorSelector* icon_color = dlg.findChild<ColorSelector*>("icon_color");
		ColorSelector* bg_color = dlg.findChild<ColorSelector*>("bg_color");
		ColorSelector* selected_color = dlg.findChild<ColorSelector*>("selected_color");
		
		peak_color->setColor(QColor(getCurrentLayer_().param.getValue("peak_color").toQString()));
		icon_color->setColor(QColor(getCurrentLayer_().param.getValue("icon_color").toQString()));
		bg_color->setColor(QColor(param_.getValue("background_color").toQString()));
		selected_color->setColor(QColor(param_.getValue("highlighted_peak_color").toQString()));
		
		if (dlg.exec())
		{
			getCurrentLayer_().param.setValue("peak_color",peak_color->getColor().name().toAscii().data());
			getCurrentLayer_().param.setValue("icon_color",icon_color->getColor().name().toAscii().data());
			param_.setValue("background_color",bg_color->getColor().name().toAscii().data());
			param_.setValue("highlighted_peak_color",selected_color->getColor().name().toAscii().data());
			
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);		
		}
	}


}//Namespace




