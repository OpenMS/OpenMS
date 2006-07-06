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
// $Id: Spectrum3DCanvas.C,v 1.25 2006/06/09 11:47:50 cfriedle Exp $
// $Author: cfriedle $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

//STL
#include <math.h>
#include<iostream.h>

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum3DCanvasPDP.h>

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
	currentDataSet().sortSpectra(true);
	currentDataSet().updateRanges(1);	
	recalculateRanges_(1,0,2);
	visible_area_.assign(overall_data_range_);
	disp_ints_.push_back(pair<float,float>(low_intensity_cutoff, overall_data_range_.max_[2]));
	emit layerActivated(this);
	openglwidget()->recalculateDotGradient_();
	recalculate_ = true;
	invalidate_();
	return current_data_;
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
		openglwidget()->setZoomFactor(1.25);
		openglwidget()->resetAngels();
		invalidate_();
		break;
	case AM_ZOOM:
		openglwidget()->setZoomFactor(1.25);
		openglwidget()->resetTranslation();
		invalidate_();
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
	openglwidget()->initializeGL();
 	openglwidget()->updateGL();
}

void Spectrum3DCanvas::intensityModeChange_()
{
	openglwidget()->recalculateDotGradient_();
	recalculate_ = true;
	invalidate_();
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
	invalidate_();
}

Spectrum3DOpenGLCanvas* Spectrum3DCanvas::openglwidget()
{
	return static_cast<Spectrum3DOpenGLCanvas*>(openglcanvas_);
}


////preferences////////////////////
  
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

