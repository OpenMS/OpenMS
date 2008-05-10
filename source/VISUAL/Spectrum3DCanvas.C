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

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum3DPrefDialog.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/FileWatcher.h>

#include <QtGui/QResizeEvent>
#include <QtGui/QComboBox>
#include <QtGui/QSpinBox>
#include <QtGui/QMenu>
#include <QtGui/QFileDialog>
#include <QtGui/QMessageBox>

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	
	Spectrum3DCanvas::Spectrum3DCanvas(const Param& preferences, QWidget* parent)
		: SpectrumCanvas(preferences, parent)
	{  
    //Paramater handling
    defaults_.setValue("dot:shade_mode", 1,"Shade mode: single-color ('flat') or gradient peaks ('smooth').");
    defaults_.setMinInt("dot:shade_mode",0);
    defaults_.setMaxInt("dot:shade_mode",1);
    defaults_.setValue("dot:gradient", "Linear|0,#ffea00;6,#ff0000;14,#aa00ff;23,#5500ff;100,#000000", "Peak color gradient.");
    defaults_.setValue("dot:interpolation_steps",200, "Interpolation steps for peak color gradient precalculation.");
    defaults_.setMinInt("dot:interpolation_steps",1);
    defaults_.setMaxInt("dot:interpolation_steps",1000);
    defaults_.setValue("dot:line_width",2,"Line width for peaks.");
    defaults_.setMinInt("dot:line_width",1);
    defaults_.setMaxInt("dot:line_width",99);
    defaults_.setValue("background_color", "#ffffff","Background color");
		setName("Spectrum3DCanvas");
		defaultsToParam_();
		setParameters(preferences);

		setFocusPolicy(Qt::TabFocus);
		openglcanvas_= new Spectrum3DOpenGLCanvas(this, *this);
		setFocusProxy(openglcanvas_);
		connect(this,SIGNAL(actionModeChange()),openglcanvas_,SLOT(actionModeChange()));
		legend_shown_ = true;
	}
		
	Spectrum3DCanvas::~Spectrum3DCanvas()
	{
	
	}
	
	void Spectrum3DCanvas::resizeEvent(QResizeEvent *e)
	{
		openglcanvas_->resize(e->size().width(),e->size().height());
		openglcanvas_->initializeGL();
	}
	
	void Spectrum3DCanvas::showLegend(bool show)
	{
		legend_shown_ = show;
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}

	bool Spectrum3DCanvas::isLegendShown() const
	{
		return legend_shown_;
	}
	
	Int Spectrum3DCanvas::finishAdding()
	{
		if (layers_.back().type!=LayerData::DT_PEAK)
		{
			QMessageBox::critical(this,"Error","This widget supports peak data only. Aborting!");
			return -1;
		}
		
		current_layer_ = getLayerCount()-1;
		currentPeakData_().sortSpectra(true);
		currentPeakData_().updateRanges(1);	

		//Abort if no data points are contained
		if (getCurrentLayer().peaks.size()==0 || getCurrentLayer().peaks.getSize()==0)
		{
			layers_.resize(getLayerCount()-1);
			if (current_layer_!=0) current_layer_ = current_layer_-1;
			QMessageBox::critical(this,"Error","Cannot add empty dataset. Aborting!");
			return -1;
		}
		
		recalculateRanges_(1,0,2);
	
		resetZoom(false);
		
		emit layerActivated(this);
		openglwidget()->recalculateDotGradient_(current_layer_);
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
		openglwidget()->updateGL();
		openglwidget()->initializeGL();

		//set watch on the file
		if (File::exists(getCurrentLayer().filename))
		{
			watcher_->addPath(getCurrentLayer().filename.toQString());
		}

		return current_layer_;
	}

	void Spectrum3DCanvas::activateLayer(int layer_index)
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()) || layer_index==int(current_layer_))
		{
			return ;
		}
		current_layer_ = layer_index;
		emit layerActivated(this);
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum3DCanvas::removeLayer(int layer_index)
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()))
		{
			return;
		}
		
		layers_.erase(layers_.begin()+layer_index);
		
		//update current layer if it became invalid
		if (current_layer_!=0 && current_layer_ >= getLayerCount()) current_layer_ = getLayerCount()-1;
		
		recalculateRanges_(1,0,2);
		
		if (layers_.empty()) return;
				
		resetZoom();
	}
	
	Spectrum3DOpenGLCanvas* Spectrum3DCanvas::openglwidget()
	{
		return static_cast<Spectrum3DOpenGLCanvas*>(openglcanvas_);
	}
	
	void Spectrum3DCanvas::update_(const char*
#ifdef DEBUG_UPDATE_
			caller_name)
	{
		cout << "Spectrum3DCanvas::update_ from '" << caller_name << "'" << endl;
#else
		)
	{
#endif
		openglwidget()->resizeGL(width(), height());
		if(update_buffer_)
		{
			update_buffer_ = false;
			if(intensity_mode_ == SpectrumCanvas::IM_SNAP)
			{
				openglwidget()->updateIntensityScale();
			}
			openglwidget()->initializeGL(); 
		}
		openglwidget()->resizeGL(width(), height());
		openglwidget()->glDraw();
	}

	void Spectrum3DCanvas::showCurrentLayerPreferences()
	{
		Internal::Spectrum3DPrefDialog dlg(this);

//		cout << "IN: " << param_ << endl;

		ColorSelector* bg_color = dlg.findChild<ColorSelector*>("bg_color");
		QComboBox* shade = dlg.findChild<QComboBox*>("shade");
		MultiGradientSelector* gradient = dlg.findChild<MultiGradientSelector*>("gradient");
		QSpinBox* width  = dlg.findChild<QSpinBox*>("width");
		QComboBox* on_file_change = dlg.findChild<QComboBox*>("on_file_change");
		
		bg_color->setColor(QColor(param_.getValue("background_color").toQString()));		
		shade->setCurrentIndex(getCurrentLayer().param.getValue("dot:shade_mode"));
		gradient->gradient().fromString(getCurrentLayer().param.getValue("dot:gradient"));
		width->setValue(UInt(getCurrentLayer().param.getValue("dot:line_width")));
		on_file_change->setCurrentIndex(on_file_change->findText(param_.getValue("on_file_change").toQString()));	

		if (dlg.exec())
		{
			param_.setValue("background_color",bg_color->getColor().name().toAscii().data());
			getCurrentLayer_().param.setValue("dot:shade_mode",shade->currentIndex());
			getCurrentLayer_().param.setValue("dot:gradient",gradient->gradient().toString());
			getCurrentLayer_().param.setValue("dot:line_width",width->value());
			param_.setValue("on_file_change", on_file_change->currentText().toAscii().data());
			
			currentLayerParamtersChanged_();
		}
	}
	
	void Spectrum3DCanvas::currentLayerParamtersChanged_()
	{
		openglwidget()->recalculateDotGradient_(current_layer_);
	 	recalculateRanges_(1,0,2);
	 	
	 	update_buffer_ = true;	
		update_(__PRETTY_FUNCTION__);
	}

	void Spectrum3DCanvas::contextMenuEvent(QContextMenuEvent* e)
	{
		QMenu* context_menu = new QMenu(this);
		QAction* result = 0;

		//Display name and warn if current layer invisible
		String layer_name = String("Layer: ") + getCurrentLayer().name;
		if (!getCurrentLayer().visible)
		{
			layer_name += " (invisible)";
		}
		context_menu->addAction(layer_name.toQString());
		context_menu->addSeparator();

		QMenu* settings_menu = new QMenu("Settings");
		settings_menu->addAction("Show/hide grid lines");
		settings_menu->addAction("Show/hide axis legends");
		settings_menu->addAction("Preferences");

		QMenu* save_menu = new QMenu("Save");
		save_menu->addAction("Layer");
		save_menu->addAction("Visible layer data");
		
		context_menu->addMenu(save_menu);
		context_menu->addMenu(settings_menu);

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
		}		
		e->accept();
	}

	void Spectrum3DCanvas::saveCurrentLayer(bool visible)
	{
  	QString file_name = QFileDialog::getSaveFileName(this, "Save file", param_.getValue("default_path").toQString(),"mzData files (*.mzData);;All files (*.*)");
		if (!file_name.isEmpty())
		{
	  	if (visible) //only visible data
	  	{
				const LayerData& layer = getCurrentLayer();
		  	DoubleReal min_mz = getVisibleArea().min()[1];
		  	DoubleReal max_mz = getVisibleArea().max()[1];
	
    		//Extract selected visible data to out
    		LayerData::ExperimentType out;
    		out.ExperimentalSettings::operator=(layer.peaks);
    		LayerData::ExperimentType::ConstIterator begin = layer.peaks.RTBegin(getVisibleArea().min()[0]);
    		LayerData::ExperimentType::ConstIterator end = layer.peaks.RTEnd(getVisibleArea().max()[0]); 
    		out.resize(end-begin);
				
				UInt i = 0;
    		for (LayerData::ExperimentType::ConstIterator it=begin; it!=end; ++it)
    		{
  				out[i].SpectrumSettings::operator=(*it);
  				out[i].setRT(it->getRT());
  				out[i].setMSLevel(it->getMSLevel());
  				out[i].setPrecursorPeak(it->getPrecursorPeak());
  				for (LayerData::ExperimentType::SpectrumType::ConstIterator it2 = it->MZBegin(min_mz); it2!= it->MZEnd(max_mz); ++it2)
  				{
  					if (layer.filters.passes(*it2))
  					{
  						out[i].push_back(*it2);
  					}
  				}
  				++i;
    		}
			  MzDataFile f;
			  f.setLogType(ProgressLogger::GUI);
			  f.store(file_name.toAscii().data(),out);
			}
			else //all data
			{
				MzDataFile().store(file_name.toAscii().data(),getCurrentLayer().peaks);
			}
		}
	}


	void Spectrum3DCanvas::updateLayer_(UInt i)
	{
		//TODO Empty layer, invalid file
		LayerData& layer = getLayer_(i);
		
		//update data
		FileHandler().loadExperiment(layer.filename,layer.peaks);
		layer.peaks.sortSpectra(true);
		layer.peaks.updateRanges(1);
		
		recalculateRanges_(1,0,2);
		resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway
		
		openglwidget()->recalculateDotGradient_(i);
		
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
		openglwidget()->updateGL();
		openglwidget()->initializeGL();
	}

	void Spectrum3DCanvas::translateLeft_()
	{
		openglwidget()->trans_x_ -= 10;
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum3DCanvas::translateRight_()
	{
		openglwidget()->trans_x_ += 10;
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum3DCanvas::translateForward_()
	{
		openglwidget()->trans_y_ += 10;
		update_(__PRETTY_FUNCTION__);		
	}
	
	void Spectrum3DCanvas::translateBackward_()
	{
		openglwidget()->trans_y_ -= 10;
		update_(__PRETTY_FUNCTION__);		
	}

}//namspace

