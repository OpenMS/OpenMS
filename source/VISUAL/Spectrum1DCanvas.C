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
// $Maintainer: Marc Sturm $
// $Authors: $
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
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DPrefDialog.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/SYSTEM/FileWatcher.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/Annotation1DItem.h>
#include <OpenMS/VISUAL/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/Annotation1DTextItem.h>
#include <OpenMS/VISUAL/Annotation1DPeakItem.h>
#include <OpenMS/VISUAL/Annotations1DContainer.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignment.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectrumAlignmentScore.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	using namespace Internal;
		
	Spectrum1DCanvas::Spectrum1DCanvas(const Param& preferences, QWidget* parent)
		: SpectrumCanvas(preferences, parent),
			mirror_mode_(false),
			show_alignment_(false),
			moving_annotations_(false)
	{
    //Paramater handling
    defaults_.setValue("highlighted_peak_color", "#ff0000", "Highlighted peak color.");
    defaults_.setValue("icon_color", "#000000", "Peak icon color.");
    defaults_.setValue("peak_color", "#0000ff", "Peak color.");
   	defaults_.setValue("annotation_color", "#000055", "Annotation color.");
    defaults_.setValue("background_color", "#ffffff", "Background color.");
		defaultsToParam_();
		setName("Spectrum1DCanvas");
		setParameters(preferences);
		
		//connect preferences change to the right slot
		connect(this,SIGNAL(preferencesChange()),this,SLOT(currentLayerParamtersChanged_()));
	}

	Spectrum1DCanvas::~Spectrum1DCanvas()
	{
	}

	void Spectrum1DCanvas::activateLayer(Size layer_index)
	{
		if (layer_index >= getLayerCount() || layer_index==current_layer_)
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
	
	void Spectrum1DCanvas::dataToWidget(const PeakType& peak, QPoint& point, bool flipped, bool percentage)
	{
		dataToWidget(peak.getMZ(), peak.getIntensity(), point, flipped, percentage);
	}
	
	void Spectrum1DCanvas::dataToWidget(float x, float y, QPoint& point, bool flipped, bool percentage)
	{
		QPoint tmp;
		if (percentage)
		{
			y *= getSnapFactor()*percentage_factor_;
		}
		SpectrumCanvas::dataToWidget_(x, y, tmp);
		point.setX(tmp.x());
		DoubleReal alignment_shrink_factor = 1.0;
		if (height() > 10)
		{
			alignment_shrink_factor = (DoubleReal)(height() - 10) / (DoubleReal)height();
		}
		if (mirror_mode_)
		{
			if (flipped)
			{
				if (!show_alignment_)
				{
					point.setY(height() - (int)(tmp.y() / 2.0));
				}
				else // show_alignment_
				{
					point.setY(height() - (int)((tmp.y() * alignment_shrink_factor) / 2.0));
				}
			}
			else // !flipped
			{
				if (!show_alignment_)
				{
					point.setY((int)(tmp.y() / 2.0));
				}
				else // show_alignment_
				{
					point.setY((int)((tmp.y() * alignment_shrink_factor) / 2.0));
				}
			}
		}
		else // !mirror_mode_
		{
			point.setY((int)(tmp.y()));
		}
	}
	
	SpectrumCanvas::PointType Spectrum1DCanvas::widgetToData(const QPoint& pos, bool percentage)
	{
		return widgetToData(pos.x(), pos.y(), percentage);
	}
	
	SpectrumCanvas::PointType Spectrum1DCanvas::widgetToData(float x, float y, bool percentage)
	{
		float actual_y;
		DoubleReal alignment_shrink_factor = 1.0;
		if (height() > 10)
		{
			alignment_shrink_factor = (DoubleReal)(height() - 10) / (DoubleReal)height();
		}
		
		if (mirror_mode_)
		{
			if (y > height() / 2)
			{
				if (!show_alignment_)
				{
					actual_y = (height() - y) * 2;
				}
				else
				{
					actual_y = (height() - y) * 2 / alignment_shrink_factor;
				}
			}
			else // y <= height()/2
			{
				if (!show_alignment_)
				{
					actual_y = y * 2;
				}
				else
				{
					actual_y = y * 2 / alignment_shrink_factor;
				}
			}
		}
		else
		{
			actual_y = y;
		}
		PointType p = SpectrumCanvas::widgetToData_(x, actual_y);
		if (percentage)
		{
			p.setY(p.getY() / (getSnapFactor()*percentage_factor_));
		}
		return p;
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
		
		if (e->button() == Qt::LeftButton)
		{
			// selection/deselection of annotation items
			Annotation1DItem* item = getCurrentLayer_().getCurrentAnnotations().getItemAt(last_mouse_pos_);
			if (item)
			{
				if (!(e->modifiers() & Qt::ControlModifier))
				{
					if (!item->isSelected())
					{
						// the item becomes the only selected item
						getCurrentLayer_().getCurrentAnnotations().deselectAll();
						item->setSelected(true);
					}
					// an item was clicked -> can be moved on the canvas
					moving_annotations_ = true;
				}
				else
				{
					// ctrl pressed -> allow selection/deselection of multiple items, do not deselect others
					item->setSelected(!item->isSelected());
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
			else
			{
				// no item was under the cursor
				getCurrentLayer_().getCurrentAnnotations().deselectAll();
			}
			
			if (action_mode_ == AM_ZOOM)
			{
				rubber_band_.setGeometry(QRect(e->pos(),QSize()));
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
						if (intensity_mode_==IM_PERCENTAGE)
						{
							percentage_factor_ = overall_data_range_.max()[1]/getCurrentLayer().getCurrentSpectrum().getMaxInt();
						}
						else 
						{
							percentage_factor_ = 1.0;
						}
						dataToWidget(peak, measurement_start_point_, getCurrentLayer().flipped);
						measurement_start_point_.setY(last_mouse_pos_.y());
					}
					else
					{
						measurement_start_.clear();
					}
				}
				else // !isMzToXAxis()
				{
					if (selected_peak_.isValid())
					{
						measurement_start_ = selected_peak_;
						const ExperimentType::PeakType& peak = measurement_start_.getPeak(getCurrentLayer().peaks);
						if (intensity_mode_==IM_PERCENTAGE)
						{
							percentage_factor_ = overall_data_range_.max()[1]/getCurrentLayer().getCurrentSpectrum().getMaxInt();
						}
						else 
						{
							percentage_factor_ = 1.0;
						}
						dataToWidget(peak, measurement_start_point_, getCurrentLayer().flipped);
						measurement_start_point_.setX(last_mouse_pos_.x());
					}
					else
					{
						measurement_start_.clear();
					}
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
		PointType data_pos = widgetToData(p);				
		emit sendCursorStatus( data_pos.getX() );
		
		PeakIndex near_peak = findPeakAtPosition_(p);
		
		if(e->buttons() & Qt::LeftButton)
		{
			bool move = moving_annotations_;
			if (mirror_mode_ && (getCurrentLayer().flipped ^ (p.y() > height()/2)))
			{
				move = false;
			}
			if (move)
			{
				if (intensity_mode_==IM_PERCENTAGE)
				{
					percentage_factor_ = overall_data_range_.max()[1]/getCurrentLayer().getCurrentSpectrum().getMaxInt();
				}
				else 
				{
					percentage_factor_ = 1.0;
				}
				PointType delta = widgetToData(p, true) - widgetToData(last_mouse_pos_, true);
				
				Annotations1DContainer& ann_1d = getCurrentLayer_().getCurrentAnnotations();
				for (Annotations1DContainer::Iterator it = ann_1d.begin(); it != ann_1d.end(); ++it)
				{
					if ((*it)->isSelected())
					{
						(*it)->move(delta);
					}
				}
				update_buffer_ = true;
				update_(__PRETTY_FUNCTION__);
				last_mouse_pos_ = p;
			}
			else if (action_mode_ == AM_TRANSLATE)
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
				if (near_peak.peak != measurement_start_.peak)
				{
					selected_peak_ = near_peak;									
					last_mouse_pos_ = p;
					update_(__PRETTY_FUNCTION__);
				}
			}
			else if (action_mode_ == AM_ZOOM)
			{
				PointType pos = widgetToData(p);
				
				if (isMzToXAxis())
				{
					rubber_band_.setGeometry(QRect(last_mouse_pos_.x(), 0, p.x() - last_mouse_pos_.x(), height()).normalized());
				}
				else
				{
					rubber_band_.setGeometry(QRect(0, last_mouse_pos_.y(), width(), p.y() - last_mouse_pos_.y()).normalized());
				}
				rubber_band_.show(); //if the mouse button is pressed before the zoom key is pressed
				
				update_(__PRETTY_FUNCTION__);
			}
		}
		else if (!e->buttons()) //no buttons pressed
		{
			selected_peak_ = findPeakAtPosition_(p);
			update_(__PRETTY_FUNCTION__);
		}
		
		//show coordinates
		if (selected_peak_.isValid())
		{
				String status;
				const ExperimentType::SpectrumType& s = selected_peak_.getSpectrum(getCurrentLayer().peaks);
				for (Size m=0; m<s.getFloatDataArrays().size();++m)
				{
					if (selected_peak_.peak < s.getFloatDataArrays()[m].size())
					{
						status += s.getFloatDataArrays()[m].getName() + ": " + s.getFloatDataArrays()[m][selected_peak_.peak] + " ";
					}
				}
				for (Size m=0; m<s.getIntegerDataArrays().size();++m)
				{
					if (selected_peak_.peak < s.getIntegerDataArrays()[m].size())
					{
						status += s.getIntegerDataArrays()[m].getName() + ": " + s.getIntegerDataArrays()[m][selected_peak_.peak] + " ";
					}
				}
				for (Size m=0; m<s.getStringDataArrays().size();++m)
				{
					if (selected_peak_.peak < s.getStringDataArrays()[m].size())
					{
						status += s.getStringDataArrays()[m].getName() + ": " + s.getStringDataArrays()[m][selected_peak_.peak] + " ";
					}
				}
				emit sendStatusMessage(status, 0);
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
				if (measurement_start_.isValid() && selected_peak_.peak != measurement_start_.peak)
				{
					const ExperimentType::PeakType& peak_1 = measurement_start_.getPeak(getCurrentLayer().peaks);
					const ExperimentType::PeakType& peak_2 = selected_peak_.getPeak(getCurrentLayer().peaks);
					DoubleReal distance = peak_2.getMZ() - peak_1.getMZ();
					// add new distance item to annotations_1d of current layer
					if (intensity_mode_==IM_PERCENTAGE)
					{
						percentage_factor_ = overall_data_range_.max()[1]/getCurrentLayer().getCurrentSpectrum().getMaxInt();
					}
					else 
					{
						percentage_factor_ = 1.0;
					}
					PointType p = widgetToData(measurement_start_point_, true);
					bool peak_1_less = peak_1.getMZ() < peak_2.getMZ();
					DoubleReal start_mz = peak_1_less ? peak_1.getMZ() : peak_2.getMZ();
					DoubleReal end_mz = peak_1_less ? peak_2.getMZ() : peak_1.getMZ();
					distance = end_mz - start_mz;
					PointType start_p(start_mz, p.getY());
					PointType end_p(end_mz, p.getY());
					
					Annotation1DItem* item = new Annotation1DDistanceItem(QString::number(distance, 'f', 3), start_p, end_p);
					getCurrentLayer_().getCurrentAnnotations().push_front(item);
				}
			}
			
			ensureAnnotationsWithinDataRange_();
			moving_annotations_ = false;
			
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
			e->accept();
			getCurrentLayer_().getCurrentAnnotations().removeSelectedItems();
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
		
		// 'a' pressed && in zoom mode (ctrl pressed) => select all annotation items
		else if ((e->modifiers() & Qt::ControlModifier) && (e->key()==Qt::Key_A))
		{
			e->accept();
			getCurrentLayer_().getCurrentAnnotations().selectAll();
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
		
		else
		{
			SpectrumCanvas::keyPressEvent(e);
		}
	}

	PeakIndex Spectrum1DCanvas::findPeakAtPosition_(QPoint p)
	{
		//no layers => return invalid peak index
		if (layers_.empty()) return PeakIndex();
		
		// mirror mode and p not on same half as active layer => return invalid peak index
		if (mirror_mode_ && (getCurrentLayer().flipped ^ (p.y() > height()/2))) return PeakIndex();
		
		//reference to the current data
		SpectrumType& spectrum = getCurrentLayer_().getCurrentSpectrum();
		Size spectrum_index = getCurrentLayer_().current_spectrum;
		
		// get the interval (in diagramm metric) that will be projected on screen coordinate p.x() or p.y() (depending on orientation)
		PointType lt = widgetToData(p - QPoint(2, 2), true);
		PointType rb = widgetToData(p + QPoint(2, 2), true);
	
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
			return PeakIndex(spectrum_index,left_it-spectrum.begin());
		}
	
		SpectrumIteratorType nearest_it = left_it;
		
		// select source interval start and end depending on diagram orientation
		if (intensity_mode_==IM_PERCENTAGE)
		{
			percentage_factor_ = overall_data_range_.max()[1]/getCurrentLayer().getCurrentSpectrum().getMaxInt();
		}
		else 
		{
			percentage_factor_ = 1.0;
		}		
		QPoint tmp;
		dataToWidget(0, overall_data_range_.minY(),tmp, getCurrentLayer().flipped, true);
		double dest_interval_start = tmp.y();
		dataToWidget(0, overall_data_range_.maxY(),tmp, getCurrentLayer().flipped, true);
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
		return PeakIndex(spectrum_index,nearest_it-spectrum.begin());
	}
	
	//////////////////////////////////////////////////////////////////////////////////
	// SLOTS
	
	void Spectrum1DCanvas::removeLayer(Size layer_index)
	{
		if (layer_index >= getLayerCount())
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
		
		if (!flippedLayersExist())
		{
			setMirrorModeActive(false);
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
			AreaType new_area(overall_data_range_.minX() - 1.0, overall_data_range_.minY(),
												overall_data_range_.maxX() + 1.0, overall_data_range_.maxY());
			changeVisibleArea_(new_area, true, true);
		}
		else
		{
			AreaType new_area(overall_data_range_.minX(), overall_data_range_.minY(),
												overall_data_range_.maxX(), overall_data_range_.maxY());
			changeVisibleArea_(new_area, true, true);
		}
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}

	void Spectrum1DCanvas::setDrawMode(DrawModes mode)
	{
		//no layers
		if (layers_.size()==0) return;
			
		if (draw_modes_[current_layer_]!=mode)
		{
			draw_modes_[current_layer_] = mode;
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
		}
	}

	Spectrum1DCanvas::DrawModes Spectrum1DCanvas::getDrawMode() const
	{ 
		//no layers
		if (layers_.size()==0) return DM_PEAKS;
			
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
		
		QTime timer;
		if (show_timing_)
		{
			timer.start();
		}
		
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
			for (Size i=0; i< getLayerCount();++i)
			{
				const LayerData& layer = getLayer(i);
				const ExperimentType::SpectrumType& spectrum = layer.getCurrentSpectrum();
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
					vbegin = getLayer_(i).getCurrentSpectrum().MZBegin(visible_area_.minX());
					vend = getLayer_(i).getCurrentSpectrum().MZEnd(visible_area_.maxX());
					// draw dashed elongations for pairs of peaks annotated with a distance
					for(Annotations1DContainer::ConstIterator it = layer.getCurrentAnnotations().begin(); it != layer.getCurrentAnnotations().end(); ++it)
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
								}
							}
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
							}
							break;
						default:
							throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
					}
					
					//draw all annotation items
					drawAnnotations(getLayer_(i), painter);
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
		
		if (mirror_mode_)
		{
			painter.save();
			
			if (!show_alignment_)
			{
				// draw x-axis
				painter.setPen(Qt::black);
				painter.drawLine(0,height()/2,width(),height()/2);
			}
			else
			{
				drawAlignment(painter);
				// two x-axes:
				painter.setPen(Qt::black);
				painter.drawLine(0,height()/2+5,width(),height()/2+5);
				painter.drawLine(0,height()/2-5,width(),height()/2-5);
			}
			
			painter.restore();
		}
		
		// draw measuring line when in measure mode and valid measurement start peak selected
		if (action_mode_ == AM_MEASURE && measurement_start_.isValid())
		{
			QPoint measurement_end_point(last_mouse_pos_.x(), measurement_start_point_.y());
			painter.drawLine(measurement_start_point_, measurement_end_point);
		}
		// draw highlighted measurement start peak and selected peak
		bool with_elongation = (action_mode_ == AM_MEASURE);
		drawHighlightedPeak_(current_layer_, measurement_start_, painter, with_elongation);
		drawHighlightedPeak_(current_layer_, selected_peak_, painter, with_elongation);
		
		//draw delta for measuring
		if (action_mode_==AM_MEASURE && measurement_start_.isValid() && selected_peak_.isValid())
		{
			drawDeltas_(painter, measurement_start_, selected_peak_, false);
		}
		else
		{
			drawCoordinates_(painter, selected_peak_, false);
		}
		
		painter.end();
#ifdef DEBUG_TOPPVIEW
		cout << "END   " << __PRETTY_FUNCTION__ << endl;
#endif
		if (show_timing_)
		{
			cout << "paint event took " << timer.elapsed() << " ms" << endl;
		}
	}
	
	void Spectrum1DCanvas::drawHighlightedPeak_(Size layer_index, const PeakIndex& peak, QPainter& painter, bool draw_elongation)
	{
		if (peak.isValid())
		{
			QPoint begin;
			const ExperimentType::PeakType& sel = peak.getPeak(getLayer_(layer_index).peaks);

			painter.setPen(QPen(QColor(param_.getValue("highlighted_peak_color").toQString()), 2));
			
			if (intensity_mode_==IM_PERCENTAGE)
			{
				percentage_factor_ = overall_data_range_.max()[1]/getLayer_(layer_index).getCurrentSpectrum().getMaxInt();
			}
			else 
			{
				percentage_factor_ = 1.0;
			}
			dataToWidget(sel, begin, getLayer_(layer_index).flipped);
			QPoint top_end(begin);
			
			bool layer_flipped = getLayer_(layer_index).flipped;
			if (isMzToXAxis())
			{
				if (layer_flipped)
				{
					top_end.setY(height());
				}
				else
				{
					top_end.setY(0);
				}
			}
			else
			{
				if (!layer_flipped)
				{
					top_end.setX(width());
				}
				else // should not happen
				{
					top_end.setX(0);
				}
			}
			
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
	
	void Spectrum1DCanvas::drawAnnotations(LayerData& layer, QPainter& painter)
	{
		bool flipped = layer.flipped;
		if (intensity_mode_==IM_PERCENTAGE)
		{
			percentage_factor_ = overall_data_range_.max()[1]/layer.getCurrentSpectrum().getMaxInt();
		}
		else 
		{
			percentage_factor_ = 1.0;
		}
		QPen pen(QColor(layer.param.getValue("annotation_color").toQString()));
		QPen selected_pen;
		
		//make selected items a little brighter
		int sel_red = pen.color().red() + 50;
		int sel_green = pen.color().green() + 50;
		int sel_blue = pen.color().blue() + 50;
		//check if rgb out of bounds
		sel_red = sel_red > 255 ? 255 : sel_red;
		sel_green = sel_green > 255 ? 255 : sel_green;
		sel_blue = sel_blue > 255 ? 255 : sel_blue;
		
		selected_pen.setColor(QColor(sel_red, sel_green, sel_blue));
		
		Annotations1DContainer& c = layer.getCurrentAnnotations();
		for (Annotations1DContainer::ConstIterator it = c.begin(); it != c.end(); ++it)
		{
			if (!(*it)->isSelected())
			{
				painter.setPen(pen);
			}
			else
			{
				painter.setPen(selected_pen);
			}
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
			QMessageBox::critical(this,"Error","Cannot add a dataset that contains no survey scans. Aborting!");
			return false;
		}
		
		//add new draw mode
		draw_modes_.push_back(DM_PEAKS);
		//estimate peak type
		PeakTypeEstimator pte;
		if (pte.estimateType(getCurrentLayer_().getCurrentSpectrum().begin(),getCurrentLayer_().getCurrentSpectrum().end()) == SpectrumSettings::RAWDATA)
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
				getCurrentLayer_().param.setValue("annotation_color", "#005500");
				break;
			case 2:
				getCurrentLayer_().param.setValue("peak_color", "#ff00ff");
				getCurrentLayer_().param.setValue("annotation_color", "#550055");
				break;
			case 3:
				getCurrentLayer_().param.setValue("peak_color", "#00ffff");
				getCurrentLayer_().param.setValue("annotation_color", "#005555");
				break;
			case 4:
				getCurrentLayer_().param.setValue("peak_color", "#ffaa00");
				getCurrentLayer_().param.setValue("annotation_color", "#550000");
				break;
		}
	
		// sort spectra in accending order of position
		for (Size i = 0; i < currentPeakData_().size(); ++i)
		{
			getCurrentLayer_().peaks[i].sortByPosition();
		}
		
		getCurrentLayer_().annotations_1d.resize(currentPeakData_().size());
		
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
		
		//Warn if negative intensities are contained
		if (getMinIntensity(current_layer_)<0.0)
		{
			QMessageBox::warning(this,"Warning","This dataset contains negative intensities. Use it at your own risk!");
		}
		
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
			DoubleReal local_max  = -numeric_limits<double>::max();
			for (Size i=0; i<getLayerCount();++i)
			{
				SpectrumType& spectrum = getLayer_(i).getCurrentSpectrum();
				SpectrumIteratorType tmp  = max_element(spectrum.MZBegin(visible_area_.minX()), spectrum.MZEnd(visible_area_.maxX()), PeakType::IntensityLess());
				if (tmp != spectrum.end() && tmp->getIntensity() > local_max) 
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
			
		  emit preferencesChange();
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
		
		Annotations1DContainer& annots_1d = getCurrentLayer_().getCurrentAnnotations();
		Annotation1DItem* annot_item = annots_1d.getItemAt(e->pos());
		if (annot_item)
		{
			annots_1d.deselectAll();
			annots_1d.selectItemAt(e->pos());
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
			
			context_menu->addAction("Edit");
			context_menu->addAction("Delete");
			if ((result = context_menu->exec(mapToGlobal(e->pos()))))
			{
				if (result->text() == "Delete")
				{
					annots_1d.removeSelectedItems();
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
	
			new_action = context_menu->addAction("Add label");
			if (mirror_mode_ && (getCurrentLayer().flipped ^ (e->pos().y() > height()/2)))
			{
				new_action->setEnabled(false);
			}
			new_action = context_menu->addAction("Add peak annotation");
			PeakIndex near_peak = findPeakAtPosition_(e->pos());
			if (!near_peak.isValid())
			{
				new_action->setEnabled(false);
			}
			context_menu->addSeparator();
			new_action = context_menu->addAction("Reset alignment");
			if (!show_alignment_)
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
						if (intensity_mode_==IM_PERCENTAGE)
						{
							percentage_factor_ = overall_data_range_.max()[1]/getCurrentLayer().getCurrentSpectrum().getMaxInt();
						}
						else 
						{
							percentage_factor_ = 1.0;
						}
						PointType position = widgetToData(e->pos(), true);
						Annotation1DItem* item = new Annotation1DTextItem(position, text);
						getCurrentLayer_().getCurrentAnnotations().push_front(item);
						
						update_buffer_ = true;
						update_(__PRETTY_FUNCTION__);
					}
				}
				else if (result->text()=="Add peak annotation")
				{
					bool ok;
					QString text = QInputDialog::getText(this, "Add peak annotation", "Enter text:", QLineEdit::Normal, "", &ok);
					if (ok && !text.isEmpty())
					{
						PeakType peak = near_peak.getPeak(getCurrentLayer().peaks);
						PointType position(peak.getMZ(), peak.getIntensity());
						Annotation1DItem* item = new Annotation1DPeakItem(position, text);
						getCurrentLayer_().getCurrentAnnotations().push_front(item);
						update_buffer_ = true;
						update_(__PRETTY_FUNCTION__);
					}
				}
				else if (result->text()=="Reset alignment")
				{
					resetAlignment();
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
		
		QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(),"mzML files (*.mzML);;All files (*)");

		if (!file_name.isEmpty())
		{
			if (visible)
			{
				ExperimentType out;
				getVisiblePeakData(out);
				addDataProcessing_(out, DataProcessing::FILTERING);
				MzMLFile().store(file_name,out);
		  }
		  else
		  {
				MzMLFile().store(file_name,layer.peaks);
		  }
		}
	}
	
	bool Spectrum1DCanvas::flippedLayersExist()
	{
		bool if_this_variable_is_true_then_there_are_flipped_layers_otherwise_not = false;
		for (Size i = 0; i < getLayerCount(); ++i)
		{
			if (layers_[i].flipped)
			{
				if_this_variable_is_true_then_there_are_flipped_layers_otherwise_not = true;
				break;
			}
		}
		return if_this_variable_is_true_then_there_are_flipped_layers_otherwise_not;
	}

	void Spectrum1DCanvas::updateLayer_(Size i)
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
	
	/// Returns whether this widget is currently in mirror mode
	bool Spectrum1DCanvas::mirrorModeActive()
	{
		return mirror_mode_;
	}
	
	/// Sets whether this widget is currently in mirror mode
	void Spectrum1DCanvas::setMirrorModeActive(bool b)
	{
		mirror_mode_ = b;
		qobject_cast<Spectrum1DWidget*>(spectrum_widget_)->toggleMirrorView(b);
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum1DCanvas::paintGridLines_(QPainter& painter)
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
				if (!mirror_mode_)
				{
					painter.drawLine(xl, y, xh, y);
				}
				else
				{
					if (!show_alignment_)
					{
						painter.drawLine(xl, y/2, xh, y/2);
						painter.drawLine(xl, yl-y/2, xh, yl-y/2);
					}
					else
					{
						DoubleReal alignment_shrink_factor = 1.0;
						if (height() > 10)
						{
							alignment_shrink_factor = (DoubleReal)(height() - 10) / (DoubleReal)height();
						}
						painter.drawLine(xl,(int)((DoubleReal)(y)*alignment_shrink_factor/2.0), xh, (int)((DoubleReal)(y)*alignment_shrink_factor/2.0));
						painter.drawLine(xl,yl-(int)((DoubleReal)(y)*alignment_shrink_factor/2.0), xh, yl-(int)((DoubleReal)(y)*alignment_shrink_factor/2.0));
					}
				}
			}
		}
		
		painter.restore();
	}
	
	void Spectrum1DCanvas::performAlignment(Size layer_index_1, Size layer_index_2, const Param& param)
	{
		alignment_.clear();
		if (layer_index_1 >= getLayerCount() || layer_index_2 >= getLayerCount())
		{
			return;
		}
		LayerData& layer_1 = getLayer_(layer_index_1);
		LayerData& layer_2 = getLayer_(layer_index_2);
		const ExperimentType::SpectrumType& spectrum_1 = layer_1.getCurrentSpectrum();
		const ExperimentType::SpectrumType& spectrum_2 = layer_2.getCurrentSpectrum();
		
		SpectrumAlignment aligner;
		aligner.setParameters(param);
		
		std::vector<std::pair<Size, Size> > aligned_peaks_indices;
		aligner.getSpectrumAlignment(aligned_peaks_indices, spectrum_1, spectrum_2);

		for (Size i = 0; i < aligned_peaks_indices.size(); ++i)
		{
			DoubleReal line_begin_mz = spectrum_1[aligned_peaks_indices[i].first].getMZ();
			DoubleReal line_end_mz = spectrum_2[aligned_peaks_indices[i].second].getMZ();
			alignment_.push_back(std::make_pair(line_begin_mz, line_end_mz));
		}
		
		show_alignment_ = true;
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
		
		SpectrumAlignmentScore scorer;
		scorer.setParameters(param);
		
		alignment_score_ = scorer(spectrum_1, spectrum_2);
	}
	
	void Spectrum1DCanvas::resetAlignment()
	{
		alignment_.clear();
		qobject_cast<Spectrum1DWidget*>(spectrum_widget_)->resetAlignment();
		show_alignment_ = false;
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum1DCanvas::drawAlignment(QPainter& painter)
	{
		painter.save();
		
		//draw peak-connecting lines between the two spectra
		painter.setPen(Qt::red);
		QPoint begin_p, end_p;
		double dummy = 0.0;

		for (Size i = 0; i < getAlignmentSize(); ++i)
		{
			dataToWidget(alignment_[i].first, dummy, begin_p);
			dataToWidget(alignment_[i].second, dummy, end_p);
			painter.drawLine(begin_p.x(), height()/2-5, end_p.x(), height()/2+5);
		}
		
		painter.restore();
	}
	
	Size Spectrum1DCanvas::getAlignmentSize()
	{
		return alignment_.size();
	}
	
	DoubleReal Spectrum1DCanvas::getAlignmentScore()
	{
		return alignment_score_;
	}
	
	void Spectrum1DCanvas::intensityModeChange_()
	{
		recalculateSnapFactor_();
		ensureAnnotationsWithinDataRange_();
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum1DCanvas::ensureAnnotationsWithinDataRange_()
	{
		for (Size i = 0; i < getLayerCount(); ++i)
		{
			if (intensity_mode_==IM_PERCENTAGE)
			{
				percentage_factor_ = overall_data_range_.max()[1]/getLayer_(i).getCurrentSpectrum().getMaxInt();
			}
			else 
			{
				percentage_factor_ = 1.0;
			}
			Annotations1DContainer& ann_1d = getLayer_(i).getCurrentAnnotations();
			for (Annotations1DContainer::Iterator it = ann_1d.begin(); it != ann_1d.end(); ++it)
			{
				(*it)->ensureWithinDataRange(this);
			}
		}
	}
	
	void Spectrum1DCanvas::flipLayer(Size index)
	{
		if (index < getLayerCount())
		{
			getLayer_(index).flipped = !getLayer_(index).flipped;
		}
	}
	
	void Spectrum1DCanvas::activateSpectrum(Size index, bool repaint)
	{
		if (index < currentPeakData_().size())
		{
			getCurrentLayer_().current_spectrum = index;
			recalculateSnapFactor_();
			if (repaint)
			{
				update_buffer_ = true;
				update_(__PRETTY_FUNCTION__);
			}
		}
	}
	
}//Namespace




