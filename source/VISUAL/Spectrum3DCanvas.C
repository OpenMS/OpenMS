// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum3DPrefDialog.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/SYSTEM/FileWatcher.h>

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
    defaults_.setValue("dot:interpolation_steps",1000, "Interpolation steps for peak color gradient precalculation.");
    defaults_.setMinInt("dot:interpolation_steps",1);
    defaults_.setMaxInt("dot:interpolation_steps",1000);
    defaults_.setValue("dot:line_width",2,"Line width for peaks.");
    defaults_.setMinInt("dot:line_width",1);
    defaults_.setMaxInt("dot:line_width",99);
    defaults_.setValue("background_color", "#ffffff","Background color");
		setName("Spectrum3DCanvas");
		defaultsToParam_();
		setParameters(preferences);

		openglcanvas_= new Spectrum3DOpenGLCanvas(this, *this);
		setFocusProxy(openglcanvas_);
		connect(this,SIGNAL(actionModeChange()),openglcanvas_,SLOT(actionModeChange()));
		legend_shown_ = true;

		//connect preferences change to the right slot
		connect(this,SIGNAL(preferencesChange()),this,SLOT(currentLayerParamtersChanged_()));
	}
		
	Spectrum3DCanvas::~Spectrum3DCanvas()
	{
	}
	
	void Spectrum3DCanvas::resizeEvent(QResizeEvent *e)
	{
		openglcanvas_->resize(e->size().width(),e->size().height());
	}
	
	void Spectrum3DCanvas::showLegend(bool show)
	{
		legend_shown_ = show;		
		update_(__PRETTY_FUNCTION__);
	}

	bool Spectrum3DCanvas::isLegendShown() const
	{
		return legend_shown_;
	}
	
	bool Spectrum3DCanvas::finishAdding_()
	{
		if (layers_.back().type!=LayerData::DT_PEAK)
		{
			QMessageBox::critical(this,"Error","This widget supports peak data only. Aborting!");
			return false;
		}
		
		current_layer_ = getLayerCount()-1;

		//Abort if no data points are contained
    if (getCurrentLayer().getPeakData()->size()==0 || getCurrentLayer().getPeakData()->getSize()==0)
		{
      layers_.resize(getLayerCount()-1);
			if (current_layer_!=0) current_layer_ = current_layer_-1;
			QMessageBox::critical(this,"Error","Cannot add a dataset that contains no survey scans. Aborting!");
			return false;
		}
		
		recalculateRanges_(0,1,2);
		resetZoom(false);
		
		//Warn if negative intensities are contained
		if (getMinIntensity(current_layer_)<0.0)
		{
			QMessageBox::warning(this,"Warning","This dataset contains negative intensities. Use it at your own risk!");
		}
			
		emit layerActivated(this);
		openglwidget()->recalculateDotGradient_(current_layer_);
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);

		return true;
	}

	void Spectrum3DCanvas::activateLayer(Size layer_index)
	{
		if (layer_index >= getLayerCount() || layer_index==current_layer_)
		{
			return ;
		}
		current_layer_ = layer_index;
		emit layerActivated(this);
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum3DCanvas::removeLayer(Size layer_index)
	{
		if (layer_index >= getLayerCount())
		{
			return;
		}
		
		layers_.erase(layers_.begin()+layer_index);
		
		//update current layer if it became invalid
		if (current_layer_!=0 && current_layer_ >= getLayerCount()) current_layer_ = getLayerCount()-1;
		
		recalculateRanges_(0,1,2);
		
		if (layers_.empty())
		{
			overall_data_range_ = DRange<3>::empty;
			update_buffer_ = true;
			update_(__PRETTY_FUNCTION__);
			return;
		}
			
		resetZoom();
	}
	
	Spectrum3DOpenGLCanvas* Spectrum3DCanvas::openglwidget()
	{
		return static_cast<Spectrum3DOpenGLCanvas*>(openglcanvas_);
	}
	
  void Spectrum3DCanvas::update_(const char* caller)
	{
    #ifdef DEBUG_TOPPVIEW
      cout << "BEGIN " << __PRETTY_FUNCTION__ << " caller: " << caller << endl;
    #endif

    std::cout << caller << std::endl;
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
		LayerData& layer = getCurrentLayer_();

//		cout << "IN: " << param_ << endl;

		ColorSelector* bg_color = dlg.findChild<ColorSelector*>("bg_color");
		QComboBox* shade = dlg.findChild<QComboBox*>("shade");
		MultiGradientSelector* gradient = dlg.findChild<MultiGradientSelector*>("gradient");
		QSpinBox* width  = dlg.findChild<QSpinBox*>("width");
		
		bg_color->setColor(QColor(param_.getValue("background_color").toQString()));		
		shade->setCurrentIndex(layer.param.getValue("dot:shade_mode"));
		gradient->gradient().fromString(layer.param.getValue("dot:gradient"));
		width->setValue(UInt(layer.param.getValue("dot:line_width")));

		if (dlg.exec())
		{
			param_.setValue("background_color",bg_color->getColor().name());
			layer.param.setValue("dot:shade_mode",shade->currentIndex());
			layer.param.setValue("dot:gradient",gradient->gradient().toString());
			layer.param.setValue("dot:line_width",width->value());
			
		  emit preferencesChange();
		}
	}
	
	void Spectrum3DCanvas::currentLayerParamtersChanged_()
	{
		openglwidget()->recalculateDotGradient_(current_layer_);
	 	recalculateRanges_(0,1,2);
	 	
	 	update_buffer_ = true;	
		update_(__PRETTY_FUNCTION__);
	}

	void Spectrum3DCanvas::contextMenuEvent(QContextMenuEvent* e)
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
		context_menu->addAction("Layer meta data");

		QMenu* save_menu = new QMenu("Save");
		context_menu->addMenu(save_menu);
		save_menu->addAction("Layer");
		save_menu->addAction("Visible layer data");

		QMenu* settings_menu = new QMenu("Settings");
		context_menu->addMenu(settings_menu);
		settings_menu->addAction("Show/hide grid lines");
		settings_menu->addAction("Show/hide axis legends");
		settings_menu->addSeparator();
 		settings_menu->addAction("Preferences");

    context_menu->addAction("Switch to 2D view");

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
			else if (result->text()=="Layer meta data")
			{
				showMetaData(true);
      } else if (result->text()=="Switch to 2D view")
      {
        emit showCurrentPeaksAs2D();
      }
		}		
		e->accept();
	}

	void Spectrum3DCanvas::saveCurrentLayer(bool visible)
	{
		const LayerData& layer = getCurrentLayer();
		
		//determine proposed filename
		String proposed_name = param_.getValue("default_path");
    if (visible==false && layer.filename!="")
    {
    	proposed_name = layer.filename;
    }
    
		QString selected_filter = "";
    QString file_name = QFileDialog::getSaveFileName(this, "Save file", proposed_name.toQString(),"mzML files (*.mzML);;mzData files (*.mzData);;mzXML files (*.mzXML);;All files (*)", &selected_filter);
    if (!file_name.isEmpty())
    {
      // check whether a file type suffix has been given
      // first check mzData and mzXML then mzML
      // if the setting is at "All files"
      // mzML will be used
    	String upper_filename = file_name;
      upper_filename.toUpper();
      if (selected_filter == "mzData files (*.mzData)")
      {
        if (!upper_filename.hasSuffix(".MZDATA"))
        {
          file_name += ".mzData";
        }
      }
      else if (selected_filter == "mzXML files (*.mzXML)")
      {
        if (!upper_filename.hasSuffix(".MZXML"))
        {
          file_name += ".mzXML";
        }
      }
      else
      {
        if (!upper_filename.hasSuffix(".MZML"))
        {
          file_name += ".mzML";
        }
      }

    	if (visible) //only visible data
    	{
				ExperimentType out;
				getVisiblePeakData(out);
				addDataProcessing_(out, DataProcessing::FILTERING);
				FileHandler().storeExperiment(file_name,out,ProgressLogger::GUI);
			}
			else //all data
			{
        FileHandler().storeExperiment(file_name,*layer.getPeakData(),ProgressLogger::GUI);
			}
		}
	}

  void Spectrum3DCanvas::updateLayer(Size i)
	{
    selected_peak_.clear();
		recalculateRanges_(0,1,2);
    resetZoom(false); //no repaint as this is done in intensityModeChange_() anyway
		openglwidget()->recalculateDotGradient_(i);
    intensityModeChange_();
    modificationStatus_(i, false);
	}

	void Spectrum3DCanvas::translateLeft_()
	{
	}
	
	void Spectrum3DCanvas::translateRight_()
	{
	}
	
	void Spectrum3DCanvas::translateForward_()
	{
	}
	
	void Spectrum3DCanvas::translateBackward_()
	{
	}

}//namspace

