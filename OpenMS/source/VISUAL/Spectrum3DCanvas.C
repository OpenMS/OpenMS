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
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum3DPrefDialog.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>

#include <QtGui/QResizeEvent>
#include <QtGui/QComboBox>
#include <QtGui/QSpinBox>

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
    defaults_.setValue("dot:gradient", "Linear|0,#efef00;11,#ffaa00;32,#ff0000;55,#aa00ff;78,#5500ff;100,#000000", "Peak color gradient.");
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
		action_mode_ = AM_ZOOM;
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
			return -1;
		}
		
		current_layer_ = getLayerCount()-1;
		currentPeakData_().sortSpectra(true);
		currentPeakData_().updateRanges(1);	
		recalculateRanges_(1,0,2);
		area_ = (getCurrentLayer().peaks.getMaxRT()-getCurrentLayer().peaks.getMinRT())*(getCurrentLayer().peaks.getMaxMZ()-getCurrentLayer().peaks.getMinMZ());
	
		visible_area_.assign(overall_data_range_);
		
		emit layerActivated(this);
		openglwidget()->recalculateDotGradient_(current_layer_);
		// update_(__PRETTY_FUNCTION__);
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
		openglwidget()->updateGL();
		openglwidget()->initializeGL();
		return current_layer_;
	}
	
	void Spectrum3DCanvas::changeVisibleArea_(const AreaType& new_area, bool add_to_stack)
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
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
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
	
	void Spectrum3DCanvas::intensityModeChange_()
	{
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
	}
	
	void Spectrum3DCanvas::removeLayer(int layer_index)
	{
		if (layer_index<0 || layer_index >= int(getLayerCount()))
		{
			return;
		}
		layers_.erase(layers_.begin()+layer_index);
		
		//update current layer
		if (current_layer_!=0 && current_layer_ >= getLayerCount())
		{
		current_layer_ = getLayerCount()-1;
		}
		
		recalculateRanges_(1,0,2);
		visible_area_.assign(overall_data_range_);
		update_buffer_ = true;
		update_(__PRETTY_FUNCTION__);
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
		
		bg_color->setColor(QColor(param_.getValue("background_color").toQString()));		
		shade->setCurrentIndex(getCurrentLayer().param.getValue("dot:shade_mode"));
		gradient->gradient().fromString(getCurrentLayer().param.getValue("dot:gradient"));
		width->setValue(UInt(getCurrentLayer().param.getValue("dot:line_width")));

		if (dlg.exec())
		{
			param_.setValue("background_color",bg_color->getColor().name().toAscii().data());
			getCurrentLayer_().param.setValue("dot:shade_mode",shade->currentIndex());
			getCurrentLayer_().param.setValue("dot:gradient",gradient->gradient().toString());
			getCurrentLayer_().param.setValue("dot:line_width",width->value());
			
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
	
}//namspace

