// -*- Mode: C++; tab-width: 2; -*-
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
 
// OpenMS
#include <OpenMS/VISUAL/PeakIcon.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DPrefDialog.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/TimeStamp.h>
#include <OpenMS/SYSTEM/FileWatcher.h>

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
			if (action_mode_ == AM_ZOOM)
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
	
		if(e->buttons() & Qt::LeftButton)
		{
			if (action_mode_ == AM_TRANSLATE)
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
			else if (action_mode_ == AM_ZOOM)
			{
				PointType pos = widgetToData_(p);
	
				rubber_band_.setGeometry(last_mouse_pos_.x(), 0, p.x() - last_mouse_pos_.x(), height());
				update_(__PRETTY_FUNCTION__);
				
				emit sendCursorStatus( pos.getX() );
			}
		}
		else if (!e->buttons()) //no buttons pressed
		{
			emit sendCursorStatus();
			selected_peak_ = findPeakAtPosition_(p);
			update_(__PRETTY_FUNCTION__);
		}
	}

	
	void Spectrum1DCanvas::mouseReleaseEvent(QMouseEvent* e)
	{
		if (e->button() == Qt::LeftButton)
		{
			if (action_mode_ == AM_ZOOM)
			{
				rubber_band_.hide();
				QRect rect = rubber_band_.geometry();
				if (rect.width()!=0)
				{
					AreaType area(widgetToData_(rect.topLeft()), widgetToData_(rect.bottomRight()));
					changeVisibleArea_(area.minX(), area.maxX(), true, true);
				}
			}
		}
	}

	PeakIndex Spectrum1DCanvas::findPeakAtPosition_(QPoint p)
	{
		//no layers => return invalid peak index
		if (layers_.empty()) return PeakIndex();
		
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
			return PeakIndex();
		}
	
		if (left_it == right_it-1 )
		{
			return PeakIndex(0,left_it-spectrum.begin());
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
					
					switch (draw_modes_[i])
					{
						case DM_PEAKS:
							//-----------------------------------------DRAWING PEAKS-------------------------------------------
							
							for (SpectrumIteratorType it = vbegin; it != vend; ++it)
							{
								if (layer.filters.passes(spectrum,it-spectrum.begin()))
								{
									dataToWidget_(*it,end);
									
									SpectrumCanvas::dataToWidget_(it->getMZ(), 0.0f, begin);
									
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
									dataToWidget_(*it, begin);
						
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
									dataToWidget_(*(vbegin-1), begin);
									dataToWidget_(*(vbegin), end);
									painter.drawLine(begin, end);
								}
							
								// clipping on right side
								if (vend!=spectrum.end() && vend!=spectrum.begin())
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
		if (selected_peak_.isValid())
		{
			const ExperimentType::PeakType& sel = selected_peak_.getPeak(getCurrentLayer().peaks);

			painter.setPen(QPen(QColor(param_.getValue("highlighted_peak_color").toQString()), 2));		
			if (getDrawMode() == DM_PEAKS)
			{
				if (intensity_mode_==IM_PERCENTAGE)
				{
					percentage_factor_ = overall_data_range_.max()[1]/getCurrentLayer().peaks[0].getMaxInt();
					dataToWidget_(sel, end);
				}
				else
				{
					dataToWidget_(sel,end);
				}
				SpectrumCanvas::dataToWidget_(sel.getMZ(), 0.0f, begin);
				painter.drawLine(begin, end);
			}
			else if (getDrawMode() == DM_CONNECTEDLINES)
			{
				if (intensity_mode_==IM_PERCENTAGE)
				{
					percentage_factor_ = overall_data_range_.max()[1]/getCurrentLayer().peaks[0].getMaxInt();
					dataToWidget_(sel, begin);
				}
				else
				{
					dataToWidget_(sel, begin);
				}
				painter.drawLine(begin.x(), begin.y()-4, begin.x(), begin.y()+4);
				painter.drawLine(begin.x()-4, begin.y(), begin.x()+4, begin.y());
			}
			emit sendCursorStatus(sel.getMZ(), sel.getIntensity());
		}


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
//						dataToWidget_(*it, end);
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
	
	void Spectrum1DCanvas::changeVisibleArea_(const AreaType& new_area, bool repaint, bool add_to_stack)
	{
		//store old zoom state
		if (add_to_stack)
		{
			zoomAdd_(new_area);
		}
		
		if (new_area!=visible_area_)
		{
			visible_area_ = new_area;
			updateScrollbars_();
			recalculateSnapFactor_();
			emit visibleAreaChanged(new_area);
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
		currentPeakData_()[0].getContainer().sortByPosition();
		
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
		QComboBox* on_file_change = dlg.findChild<QComboBox*>("on_file_change");
		
		peak_color->setColor(QColor(getCurrentLayer_().param.getValue("peak_color").toQString()));
		icon_color->setColor(QColor(getCurrentLayer_().param.getValue("icon_color").toQString()));
		bg_color->setColor(QColor(param_.getValue("background_color").toQString()));
		selected_color->setColor(QColor(param_.getValue("highlighted_peak_color").toQString()));
		on_file_change->setCurrentIndex(on_file_change->findText(param_.getValue("on_file_change").toQString()));		
		
		if (dlg.exec())
		{
			getCurrentLayer_().param.setValue("peak_color",peak_color->getColor().name());
			getCurrentLayer_().param.setValue("icon_color",icon_color->getColor().name());
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
		//Abort of there are no layers
		if (layers_.empty()) return;
		
		QMenu* context_menu = new QMenu(this);
		QAction* result = 0;

		//Display name and warn if current layer invisible
		String layer_name = String("Layer: ") + getCurrentLayer().name;
		if (!getCurrentLayer().visible)
		{
			layer_name += " (invisible)";
		}
		context_menu->addAction(layer_name.toQString())->setEnabled(false);
		context_menu->addSeparator();

		context_menu->addAction("Edit metadata");

		QMenu* save_menu = new QMenu("Save");
		save_menu->addAction("Layer");
		save_menu->addAction("Visible layer data");
		
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
			else if (result->text()=="Edit metadata")
			{
				showMetaData(true);
			}
		}		
		e->accept();
	}


	void Spectrum1DCanvas::saveCurrentLayer(bool visible)
	{
		QString file_name = QFileDialog::getSaveFileName(this, "Save file", param_.getValue("default_path").toQString(),"mzData files (*.mzData);;All files (*)");

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
				MzDataFile().store(file_name,getCurrentLayer().peaks);
		  }
			
		}
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




