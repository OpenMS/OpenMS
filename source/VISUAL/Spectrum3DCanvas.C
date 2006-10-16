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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

//STL
#include <math.h>
#include<iostream.h>

//QT
#include<qimage.h>

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum3DCanvasPDP.h>
#include <OpenMS/CONCEPT/Factory.h>

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	
Spectrum3DCanvas::Spectrum3DCanvas(QWidget* parent, const char* name, WFlags f)
	: SpectrumCanvas(parent, name, f | WRepaintNoErase)
{  
	setFocusPolicy(QWidget::TabFocus);
	openglcanvas_= new Spectrum3DOpenGLCanvas(this,"openglcanvas", *this);
	connect(openglcanvas_, SIGNAL(rightButton(QPoint)), this,SLOT(showContextMenu(QPoint)) );
	action_mode_ = AM_TRANSLATE;
	legend_shown_ = true;
	show_reduced_ = false;
}
	
Spectrum3DCanvas::~Spectrum3DCanvas()
{

}

void Spectrum3DCanvas::resizeEvent(QResizeEvent *e)
{
	openglcanvas_ ->resize(e->size().width(),e->size().height());
}

void Spectrum3DCanvas::showContextMenu(QPoint p)
{
	emit contextMenu(p);
}

void Spectrum3DCanvas::showLegend(bool show)
{
	legend_shown_ = show;
	repaintAll();
}
SignedInt Spectrum3DCanvas::finishAdding(float low_intensity_cutoff)
{
	if (type_.back()==DT_FEATURE)
	{
		datasets_.resize(datasets_.size()-1);
		features_.resize(features_.size()-1);
		type_.resize(type_.size()-1);
		return -1;
	}
	
	layer_visible_.push_back(true);
	current_data_ = getDataSetCount()-1;
	currentDataSet_().sortSpectra(true);
	currentDataSet_().updateRanges(1);	
	recalculateRanges_(1,0,2);
	if(getPrefAsInt("Preferences:3D:Data:Mode")!=0)
	{
		makeReducedDataSet();
	}
	visible_area_.assign(overall_data_range_);
	disp_ints_.push_back(pair<float,float>(low_intensity_cutoff, overall_data_range_.max_[2]));
	emit layerActivated(this);
	openglwidget()->recalculateDotGradient_();
	invalidate_();
 	repaintAll();
	return current_data_;
}
void Spectrum3DCanvas::makeReducedDataSet()
{
	reduction_param_.clear();
	switch(getPrefAsInt("Preferences:3D:Data:Mode"))
	{
	case 0:
		current_data_mode_ = 0;
		show_reduced_ = false;
		recalculateRanges_(1,0,2);
		disp_ints_.clear();
		disp_ints_.push_back(pair<float,float>( overall_data_range_.min_[2], overall_data_range_.max_[2]));
	
 		break;

	case 1:
		{	
			reduction_param_.setValue("Peaksperstep", getPrefAsInt("Preferences:3D:Data:Reduction:Max"));
			show_reduced_ = true;
			reduced_datasets_.erase(reduced_datasets_.begin(),reduced_datasets_.end());
			datareducer_ = Factory<DataReducer>::create("MaxReducer");
			for(UnsignedInt i = 0; i<datasets_.size();i++)
  		{
				ExperimentType out_experiment;
				datareducer_->setParam(reduction_param_);
				datareducer_->applyReduction(datasets_[current_data_],out_experiment);
				reduced_datasets_.push_back(out_experiment);
			}
		}
		current_data_mode_ = 1;
		recalculateRanges_(1,0,2);
		disp_ints_.clear();
		disp_ints_.push_back(pair<float,float>( overall_data_range_.min_[2], overall_data_range_.max_[2]));
		break;
		
	case 2:
		{	
			reduction_param_.setValue("Rangeperstep", getPrefAsInt("Preferences:3D:Data:Reduction:Sum"));
			//		cout<<reduction_param_.getValue("Rangeperstep")<<endl;
			show_reduced_ = true;
			reduced_datasets_.erase(reduced_datasets_.begin(),reduced_datasets_.end());
			datareducer_ = Factory<DataReducer>::create("SumReducer");
			for(UnsignedInt i = 0; i<datasets_.size();i++)
				{
					ExperimentType out_experiment;
					datareducer_->setParam(reduction_param_);
					datasets_[current_data_].sortSpectra(true);
					datareducer_->applyReduction(datasets_[current_data_],out_experiment);
					reduced_datasets_.push_back(out_experiment);
				}
			
		}
		current_data_mode_ = 2;
		recalculateRanges_(1,0,2);
		disp_ints_.clear();
		disp_ints_.push_back(pair<float,float>( overall_data_range_.min_[2], overall_data_range_.max_[2]));
		break;

	}
}

void Spectrum3DCanvas::setMainPreferences(const Param& prefs)
{
	SpectrumCanvas::setMainPreferences(prefs);
	openglwidget()->gradient_.fromString(getPrefAsString("Preferences:3D:Dot:Gradient"));
}

PreferencesDialogPage * Spectrum3DCanvas::createPreferences(QWidget* parent)
{
	return new Spectrum3DCanvasPDP(this,parent);
}

void Spectrum3DCanvas::actionModeChange_()
{
	switch(action_mode_)
	{
	case AM_TRANSLATE:
		openglwidget()->setAngels(220,220,0);
		openglwidget()->setZoomFactor(1.5,false);
		repaintAll();
		break;
	case AM_ZOOM:
		openglwidget()->setAngels(1440,0,0);
		openglwidget()->resetTranslation();
		openglwidget()->setZoomFactor(1.25,false);
		repaintAll();
		break;
	case AM_MEASURE:
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		break;
	case AM_SELECT:
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		break;
	}
}

void Spectrum3DCanvas::activateDataSet(int data_set)
{
	if ((data_set >= int(getDataSetCount())) || data_set==int(current_data_))
	{
		return ;
	}
	current_data_ = data_set;
	emit layerActivated(this);
	invalidate_();
}

void Spectrum3DCanvas::invalidate_()
{				
	if(recalculate_)
	{
		if(intensity_mode_ == SpectrumCanvas::IM_SNAP)
			{
				openglwidget()->updateIntensityScale();
			}
		openglwidget()->recalculateDotGradient_();
		openglwidget()->glInit (); 
	}	
	openglwidget()->resizeGL(width(),height());
	openglwidget()->glDraw (); 
}

void Spectrum3DCanvas::repaintAll()
{
	recalculate_ = true;
	invalidate_();

}

void Spectrum3DCanvas::intensityModeChange_()
{
	repaintAll();
}

void Spectrum3DCanvas::removeDataSet(int data_set)
{
	if (data_set >= int(getDataSetCount()))
	{
		return;
	}
	datasets_.erase(datasets_.begin()+data_set);
	features_.erase(features_.begin()+data_set);
	type_.erase(type_.begin()+data_set);
	layer_visible_.erase(layer_visible_.begin()+data_set);
	disp_ints_.erase(disp_ints_.begin()+data_set);
	recalculateRanges_(1,0,2);
	visible_area_.assign(overall_data_range_);
	repaintAll();
}

Spectrum3DOpenGLCanvas* Spectrum3DCanvas::openglwidget()
{
	return static_cast<Spectrum3DOpenGLCanvas*>(openglcanvas_);
}


////preferences////////////////////
  
SignedInt Spectrum3DCanvas::getDataMode()
{
	if(prefs_.getValue("Preferences:3D:Data:Mode").isEmpty())
	{
		return 0;
	}
	return SignedInt(prefs_.getValue("Preferences:3D:Data:Mode"));
}

void Spectrum3DCanvas::setDataMode(int data_mode)
{
	if(current_data_mode_ != data_mode || 
		(reduction_param_.getValue("Peaksperstep")!=getPref("Preferences:3D:Data:Reduction:Max") && data_mode == 1) ||
		(reduction_param_.getValue("Rangeperstep")!=getPref("Preferences:3D:Data:Reduction:Sum") && data_mode == 2))
	{
		makeReducedDataSet();
		if(zoom_stack_.empty())
		{
			resetZoom();
		}
		repaintAll();
	}
}


SignedInt Spectrum3DCanvas::getDotMode()
{
	if (prefs_.getValue("Preferences:3D:Dot:Mode").isEmpty())
	{
		return 0;
	}
	return SignedInt(prefs_.getValue("Preferences:3D:Dot:Mode"));
}

void Spectrum3DCanvas::setDotGradient(const std::string& gradient)
{
	openglcanvas_->setDotGradient(gradient);
}

SignedInt Spectrum3DCanvas::getShadeMode()
{
	if(prefs_.getValue("Preferences:3D:Shade:Mode").isEmpty())
	{
		return 0;
	}
	return SignedInt(prefs_.getValue("Preferences:3D:Shade:Mode"));
}

UnsignedInt Spectrum3DCanvas::getDotInterpolationSteps()
{
	if(prefs_.getValue("Preferences:3D:InterpolationSteps").isEmpty())
	{
		return 0;
	}
	return UnsignedInt(prefs_.getValue("Preferences:3D:InterpolationSteps"));
}

}//namspace

