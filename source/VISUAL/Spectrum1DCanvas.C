// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#include <QtGui/QMenu>
#include <QtGui/QComboBox>
#include <QtGui/QFileDialog>
#include <QtGui/QInputDialog>
 
// OpenMS
#include <OpenMS/VISUAL/PeakIcon.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DPrefDialog.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/TimeStamp.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Annotation1DItem.h>
#include <OpenMS/VISUAL/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/Annotations1DContainer.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	using namespace Internal;
		
	Spectrum1DCanvas::Spectrum1DCanvas(const Param& preferences, QWidget* parent)
		: SpectrumCanvas(preferences, parent),
			mirror_mode_(false),
			lower_half_active_(true)
	{
    //Paramater handling
    defaults_.setValue("highlighted_peak_color", "#ff0000", "Highlighted peak color.");
    defaults_.setValue("icon_color", "#000000", "Peak icon color.");
    defaults_.setValue("peak_color", "#0000ff", "Peak color.");
   	defaults_.setValue("annotation_color", "#0e6805", "Annotation color.");
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
		selected_peak_.clear();
		
		emit layerActivated(this);
	}
	
	void Spectrum1DCanvas::setVisibleArea(DRange<2> range)
	{	
		changeVisibleArea_(AreaType(range.minX(), visible_area_.minY(), range.maxX(), visible_area_.maxY()));
	}
	
	void Spectrum1DCanvas::changeVisibleArea_(double lo, double hi, bool repaint, bool add_to_stack)
	{
		changeVisibleArea_(AreaType(lo, visible_area_.minY(), hi, visible_area_.maxY()), repaint, add_to_stack);
	}
	
	void Spectrum1DCanvas::dataToWidget(const PeakType& peak, QPoint& point, bool flipped)
	{
		dataToWidget(peak.getMZ(), getSnapFactor()*percentage_factor_*peak.getIntensity(), point, flipped);
	}
	
	void Spectrum1DCanvas::dataToWidget(float x, float y, QPoint& point, bool flipped)
	{
		QPoint tmp;
		SpectrumCanvas::dataToWidget_(x, y, tmp);
		if (mirror_mode_)
		{
			if (flipped)
			{
				point.setX(tmp.x());
				point.setY(height() - (tmp.y() / 2));
			}
			else
			{
				point.setX(tmp.x());
				point.setY(tmp.y() / 2);
			}
		}
		else
		{
			point.setX(tmp.x());
			point.setY(tmp.y());
		}
	}
	
	SpectrumCanvas::PointType Spectrum1DCanvas::widgetToData(const QPoint& pos)
	{
		return widgetToData(pos.x(), pos.y());
	}
	
	SpectrumCanvas::PointType Spectrum1DCanvas::widgetToData(float x, float y)
	{
		float actual_y;
		
		if (mirror_mode_)
		{
			if (y > height() / 2)
			{
				actual_y = (height() - y) * 2;
			}
			else
			{
				actual_y = y * 2;
			}
		}
		else
		{
			actual_y = y;
		}
		
		return SpectrumCanvas::widgetToData_(x, actual_y);
	}
	
	//////////////////////////////////////////////////////////////////////////////////
	// Qt events
	
	void Spectrum1DCanvas::mousePressEvent(QMouseEvent* e)
	{
		if (current_layer_ >= getLayerCount())
		{
			return;
		}
		
		// get mouse position in widget coordinates
		last_mouse_pos_ = e->pos();
	
		if (mirror_mode_)
		{
			lower_half_active_ = (last_mouse_pos_.y() > height() / 2);
		}
		
		if (e->button() == Qt::LeftButton)
		{
			if (action_mode_ == AM_ZOOM)
			{
				rubber_band_.setGeometry(e->pos().x(),e->pos().y(),0,0);
				rubber_band_.show();
			}
			else if (action_mode_ == AM_MEASURE)
			{
				if (isMzToXAxis())
				{
					if (selected_peak_.isValid())
					{
						measurement_start_ = selected_peak_;
						const ExperimentType::PeakType& peak = measurement_start_.getPeak(getCurrentLayer().peaks);
						dataToWidget(peak, measurement_start_point_, getCurrentLayer().flipped);
						measurement_start_point_.setY(last_mouse_pos_.y());
					}
					else
					{
						measurement_start_.clear();
					}
				}
				else
				{
					QMessageBox::information(this,"Not supported","Measuring is not yet supported for rotated spectra.");
				}
			}
			/* if ctrl is pressed, allow selection / deselection of multiple annotation items.
				 else,deselect all annotation items before selecting the one under the cursor (if existent). */
			if (!(e->modifiers() & Qt::ControlModifier))
			{
				getCurrentLayer().annotations_1d.deselectAll();
			}
			Annotation1DItem* item = getCurrentLayer().annotations_1d.getItemAt(last_mouse_pos_);
			if (item)
			{
				if (item->isSelected())
				{
					getCurrentLayer().annotations_1d.deselectItemAt(last_mouse_pos_);
				}
				else
				{
					getCurrentLayer().annotations_1d.selectItemAt(last_mouse_pos_);
				}
				// if item is a distance item: show distance of selected item in status bar
				Annotation1DDistanceItem* distance_item = dynamic_cast<Annotation1DDistanceItem*>(item);
				if (distance_item)
				{
					const DoubleReal start_p = distance_item->getStartPoint().getX();
					const DoubleReal end_p = distance_item->getEndPoint().getX();
					emit sendStatusMessage(QString("Measured: dMZ = %1").arg(end_p - start_p).toStdString(), 0);
				}
			}
		}
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}

	void Spectrum1DCanvas::mouseMoveEvent(QMouseEvent* e)
	{
		if (current_layer_ >= getLayerCount())
		{
			return;
		}
		
		// mouse position relative to the diagram widget
		QPoint p = e->pos();
		
		PeakIndex near_peak = findPeakAtPosition_(p);
	
		if(e->buttons() & Qt::LeftButton)
		{
			if (action_mode_ == AM_TRANSLATE)
			{
				// translation in data metric
				double shift = widgetToData(last_mouse_pos_).getX() - widgetToData(p).getX();
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
			else if (action_mode_ == AM_MEASURE)
			{
				// measure only positive distances:
				if (near_peak.peak > measurement_start_.peak)
				{
					selected_peak_ = near_peak;
									
					if (measurement_start_.isValid() && selected_peak_.isValid())
					{
						const ExperimentType::PeakType& peak_1 = measurement_start_.getPeak(getCurrentLayer().peaks);
						const ExperimentType::PeakType& peak_2 = selected_peak_.getPeak(getCurrentLayer().peaks);
						emit sendCursorStatus(peak_2.getMZ(), peak_2.getIntensity());
						emit sendStatusMessage(QString("Measured: dMZ = %1, Intensity ratio = %2").arg(peak_2.getMZ()-peak_1.getMZ()).arg(peak_2.getIntensity()/peak_1.getIntensity()).toStdString(), 0);
					}
					else
					{
						emit sendCursorStatus();
					}
					
					last_mouse_pos_ = p;
					update_(__PRETTY_FUNCTION__);
				}
			}
			else if (action_mode_ == AM_ZOOM)
			{
				PointType pos = widgetToData(p);

				if (isMzToXAxis())
				{
					rubber_band_.setGeometry(last_mouse_pos_.x(), 0, p.x() - last_mouse_pos_.x(), height());
				}
				else
				{
					rubber_band_.setGeometry(0, last_mouse_pos_.y(), width(), p.y() - last_mouse_pos_.y());
				}
				
				update_(__PRETTY_FUNCTION__);
				
				emit sendCursorStatus( pos.getX() );
			}
		}
		else if (!e->buttons()) //no buttons pressed
		{
			selected_peak_ = findPeakAtPosition_(p);
			update_(__PRETTY_FUNCTION__);
			if (selected_peak_.isValid())
			{
				const ExperimentType::PeakType& sel = near_peak.getPeak(getCurrentLayer().peaks);
				emit sendCursorStatus(sel.getMZ(), sel.getIntensity());
			}
			else
			{
				emit sendCursorStatus();
			}
		}
	}

	
	void Spectrum1DCanvas::mouseReleaseEvent(QMouseEvent* e)
	{
		if (current_layer_ >= getLayerCount())
		{
			return;
		}
		
		if (e->button() == Qt::LeftButton)
		{
			if (action_mode_ == AM_ZOOM)
			{
				rubber_band_.hide();
				QRect rect = rubber_band_.geometry();
				if (rect.width()!=0)
				{
					AreaType area(widgetToData(rect.topLeft()), widgetToData(rect.bottomRight()));
					changeVisibleArea_(area.minX(), area.maxX(), true, true);
				}
			}
			else if (action_mode_ == AM_MEASURE)
			{
				if (!selected_peak_.isValid())
				{
					measurement_start_.clear();
				}
				if (measurement_start_.isValid() && selected_peak_.peak > measurement_start_.peak)
				{
					const ExperimentType::PeakType& peak_1 = measurement_start_.getPeak(getCurrentLayer().peaks);
					const ExperimentType::PeakType& peak_2 = selected_peak_.getPeak(getCurrentLayer().peaks);
					emit sendCursorStatus(peak_2.getMZ(), peak_2.getIntensity());
					DoubleReal distance = peak_2.getMZ() - peak_1.getMZ();
					emit sendStatusMessage(QString("Measured: dMZ = %1, Intensity ratio = %2").arg(distance).arg(peak_2.getIntensity()/peak_1.getIntensity()).toStdString(), 0);
					// add new distance item to annotations_1d of current layer
					PointType start_p = widgetToData(measurement_start_point_);
					start_p.setX(peak_1.getMZ());
					PointType end_p(peak_2.getMZ(), start_p.getY());
					QPen pen(QColor(getCurrentLayer().param.getValue("annotation_color").toQString()), 1);
					Annotation1DItem* item = new Annotation1DDistanceItem(QString::number(distance, 'f', 3), start_p, end_p, pen);
					getCurrentLayer().annotations_1d.push_front(item);
				}
			}
			measurement_start_.clear();
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
	}

	void Spectrum1DCanvas::keyPressEvent(QKeyEvent* e)
	{
		// Delete pressed => delete selected annotations from the current layer
		if (e->key()==Qt::Key_Delete)
		{
			getCurrentLayer().annotations_1d.removeSelectedItems();
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
		
		// 'a' pressed && in zoom mode (ctrl pressed) => select all annotation items
		if ((e->modifiers() & Qt::ControlModifier) && (e->key()==Qt::Key_A))
		{
			getCurrentLayer().annotations_1d.selectAll();
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
		
		// Alt/Shift pressed => change action mode
		if (e->key()==Qt::Key_Control)
		{
			action_mode_ = AM_ZOOM;
			emit actionModeChange();
		}
		else if (e->key()==Qt::Key_Shift)
		{
			action_mode_ = AM_MEASURE;
			emit actionModeChange();
		}
		
		// CTRL+/CTRL- => Zoom stack
		if ((e->modifiers() & Qt::ControlModifier) && (e->key()==Qt::Key_Plus))
		{
			zoomForward_();
		}
		else if ((e->modifiers() & Qt::ControlModifier) && (e->key()==Qt::Key_Minus))
		{
			zoomBack_();
		}
		
		// Arrow keys => translate
		else if (e->key()==Qt::Key_Left)
		{
			translateLeft_();
		}
		else if (e->key()==Qt::Key_Right)
		{
			translateRight_();
		}
		else if (e->key()==Qt::Key_Up)
		{
			translateForward_();
		}
		else if (e->key()==Qt::Key_Down)
		{
			translateBackward_();
		}
		
		//Backspace to reset zoom
		else if (e->key()==Qt::Key_Backspace)
		{
			resetZoom();
		}
		
		e->ignore();
		
	}

	PeakIndex Spectrum1DCanvas::findPeakAtPosition_(QPoint p)
	{
		//no layers => return invalid peak index
		if (layers_.empty()) return PeakIndex();
		
		//reference to the current data
		SpectrumType& spectrum = currentPeakData_()[0];
		
		// get the interval (in diagramm metric) that will be projected on screen coordinate p.x() or p.y() (depending on orientation)
		PointType lt = widgetToData(p - QPoint(1, 1));
		PointType rb = widgetToData(p + QPoint(1, 1));
	
		// get iterator on first peak with higher position than interval_start
		PeakType temp;
		temp.setMZ(min(lt.getX(),rb.getX()));
		SpectrumIteratorType left_it = lower_bound(spectrum.begin(), spectrum.end(), temp, PeakType::PositionLess());
	
		// get iterator on first peak with higher position than interval_end
		temp.setMZ(max(lt.getX(),rb.getX()));
		SpectrumIteratorType	right_it = lower_bound(left_it, spectrum.end(), temp, PeakType::PositionLess());
	
	
		if (left_it == right_it) // both are equal => no peak falls into this interval
		{
			return PeakIndex();
		}
	
		if (left_it == right_it-1 )
		{
			return PeakIndex(0,left_it-spectrum.begin());
		}
	
		SpectrumIteratorType nearest_it = left_it;
		
		// select source interval start and end depending on diagram orientation
		QPoint tmp;
		dataToWidget(0, overall_data_range_.minY(),tmp, getCurrentLayer().flipped);
		double dest_interval_start = tmp.y();
		dataToWidget(0, overall_data_range_.maxY(),tmp, getCurrentLayer().flipped);
		double dest_interval_end = tmp.y();
		
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
	
		return PeakIndex(0,nearest_it-spectrum.begin());
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
	
		//update current layer if it became invalid
		if (current_layer_!=0 && current_layer_ >= getLayerCount()) current_layer_ = getLayerCount()-1;
		
		//update nearest peak
		selected_peak_.clear();
		
		//abort if there are no layers anymore
		if (layers_.empty())
		{
			overall_data_range_ = DRange<3>::empty;
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
			return;
		}
		
		//update range area
		recalculateRanges_(0,2,1);
		overall_data_range_.setMinY(0.0);  // minimal intensity always 0.0
		float width = overall_data_range_.width();
		overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
		overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
		overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());
		
		zoomClear_();
		if (overall_data_range_.maxX() - overall_data_range_.minX() <1.0)
		{
			changeVisibleArea_(overall_data_range_.minX() -1.0, overall_data_range_.maxX() + 1.0, true, true);
		}
		else
		{
			changeVisibleArea_(overall_data_range_.minX(), overall_data_range_.maxX(), true, true);
		}
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
		//Only fill background if no layer is present
		if (getLayerCount()==0)
		{
			QPainter painter;
			painter.begin(this);
			painter.fillRect(0,0,this->width(),this->height(),QColor(param_.getValue("background_color").toQString()));
			painter.end();
			e->accept();
			return;
		}

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
				const LayerData& layer = getLayer(i);
				const ExperimentType::SpectrumType& spectrum = layer.peaks[0];
				if (layer.visible)
				{
					QPen icon_pen = QPen(QColor(layer.param.getValue("icon_color").toQString()), 1);
					painter.setPen(QPen(QColor(layer.param.getValue("peak_color").toQString()), 1));
					if (intensity_mode_ == IM_PERCENTAGE)
					{
						percentage_factor_ = overall_data_range_.max()[1]/spectrum.getMaxInt();
					}
					else 
					{
						percentage_factor_ = 1.0;
					}
					
					vbegin = getLayer_(i).peaks[0].MZBegin(visible_area_.minX());
					vend = getLayer_(i).peaks[0].MZEnd(visible_area_.maxX());
					
					// draw dashed elongations for pairs of peaks annotated with a distance
					for(Annotations1DContainer::Iterator it = layer.annotations_1d.begin(); it != layer.annotations_1d.end(); ++it)
					{
						Annotation1DDistanceItem* distance_item = dynamic_cast<Annotation1DDistanceItem*>(*it);
						if (distance_item)
						{
								QPoint from;
								QPoint to;
								dataToWidget(distance_item->getStartPoint().getX(), 0, from, layer.flipped);
								
								dataToWidget(distance_item->getStartPoint().getX(), getVisibleArea().maxY(), to, layer.flipped);
								drawDashedLine_(from, to, painter);
								
								dataToWidget(distance_item->getEndPoint().getX(), 0, from, layer.flipped);
								
								dataToWidget(distance_item->getEndPoint().getX(), getVisibleArea().maxY(), to, layer.flipped);
								drawDashedLine_(from, to, painter);
						}
					}
					
					switch (draw_modes_[i])
					{
						case DM_PEAKS:
							//-----------------------------------------DRAWING PEAKS-------------------------------------------
							
							for (SpectrumIteratorType it = vbegin; it != vend; ++it)
							{
								if (layer.filters.passes(spectrum,it-spectrum.begin()))
								{
									dataToWidget(*it,end,layer.flipped);
									
									dataToWidget(it->getMZ(), 0.0f, begin, layer.flipped);
									
									// draw peak
									painter.drawLine(begin, end);

//									//draw icon if necessary
//									if (it->metaValueExists(4))
//									{
//										painter.save();
//										painter.setPen(icon_pen);	
//										PeakIcon::drawIcon((PeakIcon::Icon)(UInt)(it->getMetaValue(4)),painter,QRect(end.x() - 5, end.y() - 5, 10, 10));
//										painter.restore();
//									}
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
									dataToWidget(*it, begin, layer.flipped);
						
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
									
//									// draw associated icon
//									if (it->metaValueExists(4))
//									{
//										painter.save();
//										painter.setPen(icon_pen);											
//										PeakIcon::drawIcon((PeakIcon::Icon)(UInt)(it->getMetaValue(4)),painter,QRect(begin.x() - 5, begin.y() - 5, 10, 10));
//										painter.restore();
//									}
								}
								painter.drawPath(path);
									
								// clipping on left side
								if (vbegin!=spectrum.begin() && vbegin!=spectrum.end())
								{
									dataToWidget(*(vbegin-1), begin, layer.flipped);
									dataToWidget(*(vbegin), end, layer.flipped);
									painter.drawLine(begin, end);
								}
							
								// clipping on right side
								if (vend!=spectrum.end() && vend!=spectrum.begin())
								{
									dataToWidget(*(vend-1), begin, layer.flipped);
									dataToWidget(*(vend), end, layer.flipped);
									painter.drawLine(begin,end);
								}
								//-------------------------------------DRAWING CONNECTED LINES END-----------------------------------------
							}
							break;
						default:
							throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
					}
					
					//draw all annotation items
					drawAnnotations(layer, painter);
				}
			}
			painter.end();
		} //if (update_buffer)
		painter.begin(this);
		
		//draw peak data
		QVector<QRect> rects = e->region().rects();
		for (int i = 0; i < (int)rects.size(); ++i)
		{
			painter.drawImage(rects[i].topLeft(), buffer_, rects[i]);
		}
		
		// draw x-axis when in mirror mode
		if (mirror_mode_)
		{
			painter.save();
			painter.setPen(Qt::black);
			painter.drawLine(0,height()/2,width(),height()/2);
			painter.restore();
		}
		
		// draw measuring line when in measure mode and valid measurement start peak selected
		if (action_mode_ == AM_MEASURE && measurement_start_.isValid())
		{
			QPoint measurement_end_point(last_mouse_pos_.x(), measurement_start_point_.y());
			if (measurement_end_point.x() >= measurement_start_point_.x())
			{
				painter.drawLine(measurement_start_point_, measurement_end_point);
			}
		}
		// draw highlighted measurement start peak and selected peak
		bool with_elongation = (action_mode_ == AM_MEASURE);
		drawHighlightedPeak_(current_layer_, measurement_start_, painter, with_elongation);
		drawHighlightedPeak_(current_layer_, selected_peak_, painter, with_elongation);
		
//		if (draw_metainfo_)
//		{
//			SpectrumIteratorType vbegin, vend;
//			for (UInt i=0; i< getLayerCount();++i)
//			{
//				if (getLayer(i).visible)
//				{
//
//					vbegin = getLayer_(i).peaks[0].MZBegin(visible_area_.minX());
//					vend = getLayer_(i).peaks[0].MZEnd(visible_area_.maxX());
//			
//					for (SpectrumIteratorType it = vbegin; it != vend; it++)
//					{
//						dataToWidget(*it, end);
//						painter.drawText(end, it->getMetaValue("IonName").toQString());
//					}
//				}
//			}
//		}

		painter.end();
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
#ifdef TIMING_TOPPVIEW	
		cout << "1D PaintEvent took " << timer.elapsed() << " ms" << endl;
#endif	

	}
	
	void Spectrum1DCanvas::drawHighlightedPeak_(UInt layer_index, const PeakIndex& peak, QPainter& painter, bool draw_elongation)
	{
		if (peak.isValid())
		{
			QPoint begin;
			const ExperimentType::PeakType& sel = peak.getPeak(getLayer_(layer_index).peaks);

			painter.setPen(QPen(QColor(param_.getValue("highlighted_peak_color").toQString()), 2));
				
			if (intensity_mode_==IM_PERCENTAGE)
			{
				percentage_factor_ = overall_data_range_.max()[1]/getLayer_(layer_index).peaks[0].getMaxInt();
			}
			
			dataToWidget(sel, begin, getLayer_(layer_index).flipped);
			QPoint top_end;
			dataToWidget(sel.getMZ(), getVisibleArea().maxY(), top_end, getLayer_(layer_index).flipped);					
			
			// paint the crosshair only for currently selected peaks of the current layer
			if (layer_index == current_layer_ && (peak == measurement_start_ || peak == selected_peak_))
			{
				painter.drawLine(begin.x(), begin.y()-4, begin.x(), begin.y()+4);
				painter.drawLine(begin.x()-4, begin.y(), begin.x()+4, begin.y());
			}
			// draw elongation as dashed line (while in measure mode and for all existing distance annotations)
			if (draw_elongation)
			{
				drawDashedLine_(begin, top_end, painter);
			}
		}
	}
	
	void Spectrum1DCanvas::drawDashedLine_(const QPoint& from, const QPoint& to, QPainter& painter)
	{
		QPen pen;
		QVector<qreal> dashes;
		dashes << 5 << 5 << 1 << 5;
		pen.setDashPattern(dashes);
		pen.setColor(QColor(param_.getValue("highlighted_peak_color").toQString()));
		painter.save();
		painter.setPen(pen);
		painter.drawLine(from, to);
		painter.restore();
	}
	
	void Spectrum1DCanvas::drawAnnotations(const LayerData& layer, QPainter& painter)
	{
		bool flipped = layer.flipped;
		Annotations1DContainer& c = layer.annotations_1d;
		for (Annotations1DContainer::ConstIterator it = c.begin(); it != c.end(); ++it)
		{
			(*it)->draw(this, painter, flipped);
		}
	}
	
	void Spectrum1DCanvas::changeVisibleArea_(const AreaType& new_area, bool repaint, bool add_to_stack)
	{
		if (new_area!=visible_area_)
		{
			visible_area_ = new_area;
			updateScrollbars_();
			recalculateSnapFactor_();
			emit visibleAreaChanged(new_area);
		}
		
		//store old zoom state
		if (add_to_stack)
		{
			zoomAdd_(new_area);
		}
		
		if (repaint)
		{
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
	}
	
	// Destructor
	Spectrum1DCanvas::~Spectrum1DCanvas()
	{
		
	}
	
	bool Spectrum1DCanvas::finishAdding_()
	{
		if (layers_.back().type!=LayerData::DT_PEAK)
		{
			QMessageBox::critical(this,"Error","This widget supports peak data only. Aborting!");
			return false;
		}
		
		current_layer_ = getLayerCount()-1;
		currentPeakData_().updateRanges();
		
		//Abort if no data points are contained
		if (getCurrentLayer().peaks.size()==0 || getCurrentLayer().peaks.getSize()==0)
		{
			layers_.resize(getLayerCount()-1);
			if (current_layer_!=0) current_layer_ = current_layer_-1;
			QMessageBox::critical(this,"Error","Cannot add an empty dataset. Aborting!");
			return false;
		}
				
		//add new draw mode
		draw_modes_.push_back(DM_PEAKS);
		//estimate peak type
		PeakTypeEstimator pte;
		if (pte.estimateType(currentPeakData_()[0].begin(),currentPeakData_()[0].end()) == SpectrumSettings::RAWDATA)
		{
			draw_modes_.back() = DM_CONNECTEDLINES;
		}
		
		//Change peak color if this is not the first layer
		switch(current_layer_%5)
		{
			case 0:
				break;
			case 1:
				getCurrentLayer_().param.setValue("peak_color", "#00ff00");
				break;
			case 2:
				getCurrentLayer_().param.setValue("peak_color", "#ff00ff");
				break;
			case 3:
				getCurrentLayer_().param.setValue("peak_color", "#00ffff");
				break;
			case 4:
				getCurrentLayer_().param.setValue("peak_color", "#ffaa00");
				break;
		}
	
		// sort peaks in accending order of position
		currentPeakData_()[0].sortByPosition();
		
		//update nearest peak
		selected_peak_.clear();
		
		//update ranges
		recalculateRanges_(0,2,1);
		overall_data_range_.setMinY(0.0);  // minimal intensity always 0.0
		float width = overall_data_range_.width();
		overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
		overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
		overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());
		
		resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway
		
		if (getLayerCount()==2)
		{
			setIntensityMode(IM_PERCENTAGE);
		}
		intensityModeChange_();

		emit layerActivated(this);
		
		//set watch on the file
		if (File::exists(getCurrentLayer().filename))
		{
			watcher_->addFile(getCurrentLayer().filename.toQString());
		}
		
		return true;
	}

  void Spectrum1DCanvas::recalculateSnapFactor_()
  {
  	if (intensity_mode_ == IM_SNAP) 
		{
			double local_max  = -numeric_limits<double>::max();
			for (UInt i=0; i<getLayerCount();++i)
			{
				SpectrumIteratorType tmp  = max_element(getLayer_(i).peaks[0].MZBegin(visible_area_.minX()), getLayer_(i).peaks[0].MZEnd(visible_area_.maxX()), PeakType::IntensityLess());
				if (tmp->getIntensity() > local_max) 
				{
					local_max = tmp->getIntensity();
				}
			}
			snap_factors_[0] = overall_data_range_.max()[1]/local_max;
		}
		else
		{ 
			snap_factors_[0] = 1.0;
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
		ColorSelector* annotation_color = dlg.findChild<ColorSelector*>("annotation_color");
		ColorSelector* bg_color = dlg.findChild<ColorSelector*>("bg_color");
		ColorSelector* selected_color = dlg.findChild<ColorSelector*>("selected_color");
		QComboBox* on_file_change = dlg.findChild<QComboBox*>("on_file_change");
		
		peak_color->setColor(QColor(getCurrentLayer_().param.getValue("peak_color").toQString()));
		icon_color->setColor(QColor(getCurrentLayer_().param.getValue("icon_color").toQString()));
		annotation_color->setColor(QColor(getCurrentLayer_().param.getValue("annotation_color").toQString()));
		bg_color->setColor(QColor(param_.getValue("background_color").toQString()));
		selected_color->setColor(QColor(param_.getValue("highlighted_peak_color").toQString()));
		on_file_change->setCurrentIndex(on_file_change->findText(param_.getValue("on_file_change").toQString()));		
		
		if (dlg.exec())
		{
			getCurrentLayer_().param.setValue("peak_color",peak_color->getColor().name());
			getCurrentLayer_().param.setValue("icon_color",icon_color->getColor().name());
			getCurrentLayer_().param.setValue("annotation_color",annotation_color->getColor().name());
			param_.setValue("background_color",bg_color->getColor().name());
			param_.setValue("highlighted_peak_color",selected_color->getColor().name());
			param_.setValue("on_file_change", on_file_change->currentText());
			
			currentLayerParamtersChanged_();
		}
	}

	void Spectrum1DCanvas::currentLayerParamtersChanged_()
	{
		update_buffer_ = true;	
		update_(__PRETTY_FUNCTION__);
	}

	void Spectrum1DCanvas::contextMenuEvent(QContextMenuEvent* e)
	{
		//Abort if there are no layers
		if (layers_.empty()) return;
		
		QMenu* context_menu = new QMenu(this);
		QAction* result = 0;
		QAction* new_action = 0;

		Annotation1DItem* annot_item = getCurrentLayer().annotations_1d.getItemAt(e->pos());
		if (annot_item)
		{
			getCurrentLayer().annotations_1d.deselectAll();
			getCurrentLayer().annotations_1d.selectItemAt(e->pos());
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
			
			context_menu->addAction("Edit");
			context_menu->addAction("Delete");
			if ((result = context_menu->exec(mapToGlobal(e->pos()))))
			{
				if (result->text() == "Delete")
				{
					getCurrentLayer().annotations_1d.removeSelectedItems();
				}
				else if (result->text() == "Edit")
				{
					const String& old_text = annot_item->getText();
					
					bool ok;
					QString text = QInputDialog::getText(this, "Edit text", "Enter text:", QLineEdit::Normal, old_text.toQString(), &ok);
					if (ok && !text.isEmpty())
					{
						annot_item->setText(text);
					}
				}
				update_buffer_ = true;
				update_(__PRETTY_FUNCTION__);
			}
		}
		else
		{
			//Display name and warn if current layer invisible
			String layer_name = String("Layer: ") + getCurrentLayer().name;
			if (!getCurrentLayer().visible)
			{
				layer_name += " (invisible)";
			}
			context_menu->addAction(layer_name.toQString())->setEnabled(false);
			context_menu->addSeparator();
	
			context_menu->addAction("Add label");
			new_action = context_menu->addAction("Add peak annotation");
			PeakIndex near_peak = findPeakAtPosition_(e->pos());
			if (!near_peak.isValid())
			{
				new_action->setEnabled(false);
			}
			
			context_menu->addSeparator();
	
			context_menu->addAction("Layer meta data");
	
			QMenu* save_menu = new QMenu("Save");
			save_menu->addAction("Layer");
			save_menu->addAction("Visible layer data");
			save_menu->addAction("As image");
			
			QMenu* settings_menu = new QMenu("Settings");
			settings_menu->addAction("Show/hide grid lines");
			settings_menu->addAction("Show/hide axis legends");
			settings_menu->addAction("Show as raw data/peaks");
			settings_menu->addSeparator();
			settings_menu->addAction("Preferences");
			
			context_menu->addMenu(save_menu);
			context_menu->addMenu(settings_menu);
	
			//add external context menu
			if (context_add_)
			{
				context_menu->addSeparator();
				context_menu->addMenu(context_add_);
			}
	
			//evaluate menu
			if ((result = context_menu->exec(mapToGlobal(e->pos()))))
			{
				if (result->text() == "Preferences")
				{
					showCurrentLayerPreferences();
				}
				else if (result->text() == "Show/hide grid lines")
				{
					showGridLines(!gridLinesShown());
				} 
				else if (result->text() == "Show/hide axis legends")
				{
					emit changeLegendVisibility();
				}
				else if (result->text()=="Layer" || result->text()=="Visible layer data")
				{
					saveCurrentLayer(result->text()=="Visible layer data");
				}
				else if (result->text()=="As image")
				{
					spectrum_widget_->saveAsImage();
				}
				else if (result->text()=="Show as raw data/peaks")
				{
					if (getDrawMode()==DM_PEAKS)
					{
						setDrawMode(DM_CONNECTEDLINES);
					}
					else
					{
						setDrawMode(DM_PEAKS);
					}
				}
				else if (result->text()=="Layer meta data")
				{
					showMetaData(true);
				}
				else if (result->text()=="Add label")
				{
					bool ok;
					QString text = QInputDialog::getText(this, "Add label", "Enter text:", QLineEdit::Normal, "", &ok);
					if (ok && !text.isEmpty())
					{
						PointType position = widgetToData(e->pos());
						QPen pen(QColor(getCurrentLayer().param.getValue("annotation_color").toQString()), 1);
						Annotation1DItem* item = new Annotation1DTextItem(position, text, pen);
						getCurrentLayer().annotations_1d.push_front(item);
						
						update_buffer_ = true;
						update_(__PRETTY_FUNCTION__);
					}
				}
				else if  (result->text()=="Add peak annotation")
				{
					bool ok;
					QString text = QInputDialog::getText(this, "Add peak annotation", "Enter text:", QLineEdit::Normal, "", &ok);
					if (ok && !text.isEmpty())
					{
						PointType position = widgetToData(e->pos());
						PeakType peak = near_peak.getPeak(getCurrentLayer().peaks);
						position.setX(peak.getMZ());
						DoubleReal peak_height = peak.getIntensity();
						PointType actual_pos = position;
						if (position.getY() > peak_height)
						{
							actual_pos.setY(peak_height);
						}
						QPen pen(QColor(getCurrentLayer().param.getValue("annotation_color").toQString()), 1);
						Annotation1DItem* item = new Annotation1DPeakItem(actual_pos, near_peak, text, pen);
						getCurrentLayer().annotations_1d.push_front(item);
						update_buffer_ = true;
						update_(__PRETTY_FUNCTION__);
					}
				}
			}
		}
		e->accept();
	}


	void Spectrum1DCanvas::saveCurrentLayer(bool visible)
	{
		const LayerData& layer = getCurrentLayer();
		
		//determine proposed filename
		String proposed_name = param_.getValue("default_path");
    if (visible==false && layer.filename!="")
    {
    	proposed_name = layer.filename;
    }
		
		QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(),"mzData files (*.mzData);;All files (*)");

		if (!file_name.isEmpty())
		{
			if (visible)
			{
				ExperimentType out;
				getVisiblePeakData(out);
				MzDataFile().store(file_name,out);
		  }
		  else
		  {
				MzDataFile().store(file_name,layer.peaks);
		  }
		}
	}
	
	bool Spectrum1DCanvas::flippedLayersExist()
	{
		bool if_this_variable_is_true_then_there_are_flipped_layers_otherwise_not = false;
		for (UInt i = 0; i < layers_.size(); ++i)
		{
			if (layers_[i].flipped)
			{
				if_this_variable_is_true_then_there_are_flipped_layers_otherwise_not = true;
			}
		}
		return if_this_variable_is_true_then_there_are_flipped_layers_otherwise_not;
	}

	void Spectrum1DCanvas::updateLayer_(UInt i)
	{
		LayerData& layer = getLayer_(i);
		try
		{
			FileHandler().loadExperiment(layer.filename,layer.peaks);
		}
		catch(Exception::BaseException& e)
		{
			QMessageBox::critical(this,"Error",(String("Error while loading file") + layer.filename + "\nError message: " + e.what()).toQString());
			layer.peaks.clear();
		}		
		layer.peaks.resize(1);
		layer.peaks.sortSpectra();
		layer.peaks.updateRanges();
		
		//update nearest peak
		selected_peak_.clear();
		
		//update ranges
		recalculateRanges_(0,2,1);
		overall_data_range_.setMinY(0.0);  // minimal intensity always 0.0
		float width = overall_data_range_.width();
		overall_data_range_.setMinX(overall_data_range_.minX() - 0.002 * width);
		overall_data_range_.setMaxX(overall_data_range_.maxX() + 0.002 * width);
		overall_data_range_.setMaxY(overall_data_range_.maxY() + 0.002 * overall_data_range_.height());
		
		resetZoom();
		modificationStatus_(i, false);
	}

	void Spectrum1DCanvas::translateLeft_()
	{
		DoubleReal shift = 0.05 * visible_area_.width();
		DoubleReal newLo = visible_area_.minX() - shift;
		DoubleReal newHi = visible_area_.maxX() - shift;
		// check if we are falling out of bounds
		if (newLo < overall_data_range_.minX())
		{
			newLo = overall_data_range_.minX();
			newHi = newLo + visible_area_.width();
		}
		//chage data area
		changeVisibleArea_(newLo, newHi);
	}
	
	void Spectrum1DCanvas::translateRight_()
	{
		DoubleReal shift = 0.05 * visible_area_.width();
		DoubleReal newLo = visible_area_.minX() + shift;
		DoubleReal newHi = visible_area_.maxX() + shift;
		// check if we are falling out of bounds
		if (newHi > overall_data_range_.maxX())
		{
			newHi = overall_data_range_.maxX();
			newLo = newHi - visible_area_.width();
		}
		//chage data area
		changeVisibleArea_(newLo, newHi);
	}
}//Namespace




