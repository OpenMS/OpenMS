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
// $Id: Spectrum3DCanvasPDP.C,v 1.12 2006/05/30 13:48:19 cfriedle Exp $
// $Author: cfriedle $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/VISUAL/DIALOGS/Spectrum3DCanvasPDP.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>


// Qt
#include <qlayout.h>
#include <qradiobutton.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qvbuttongroup.h>
#include <qspinbox.h>
using namespace std;

namespace OpenMS
{
	namespace Internal
	{
		Spectrum3DCanvasPDP::Spectrum3DCanvasPDP(Spectrum3DCanvas* manager, QWidget* parent, const char* name, WFlags f)
				: PreferencesDialogPage(manager,parent,name,f)
		{
			help_ = "This is the preferences dialog of 3D spectrum!"
								"<br>";
			QGridLayout* grid;
			QLabel * label;
			grid = new QGridLayout(this, 2, 1);
			grid->setMargin(6);
			grid->setSpacing(4);	
			QGroupBox* box = new QGroupBox(2,Qt::Horizontal,"Dot coloring",this);
			QVButtonGroup* coloring_group = new QVButtonGroup("Color Mode:",box);
			//coloring_group->setFrameStyle(QFrame::NoFrame);
			box->addSpace(0);
			dot_mode_black_ = new QRadioButton("Black",coloring_group);
			dot_mode_gradient_ = new QRadioButton("Gradient",coloring_group);
			dot_gradient_ = new MultiGradientSelector(box);
			box->addSpace(0);
			label = new QLabel("Interpolation steps: ",box);
			dot_interpolation_steps_ = new QSpinBox(10,1000,1,box,"");
			dot_interpolation_steps_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			QVButtonGroup* shading_group = new QVButtonGroup("Shade Mode:",box);
			//shading_group->setFrameStyle(QFrame::NoFrame);
			shade_mode_flat_ = new QRadioButton("Flat",shading_group);
			shade_mode_smooth_ = new QRadioButton("Smooth",shading_group);			
			grid->addMultiCellWidget(box,0,1,0,0);

			box = new QGroupBox(2,Qt::Horizontal,"Line Width",this);
			label = new QLabel("Line Width: ",box);
			dot_line_width_ = new QSpinBox(1,10,1,box,"");
			dot_line_width_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			grid->addWidget(box,2,0);	

			box = new QGroupBox(2,Qt::Horizontal,"Colors",this);
			label = new QLabel("Background color: ",box);
			background_color_ = new ColorSelector(box);
			label = new QLabel("Axes Color: ",box);
			axes_color_ = new ColorSelector(box);
			grid->addWidget(box,0,1);	

			box = new QGroupBox(2,Qt::Horizontal,"Intensity",this);
			QVButtonGroup* intenstiy_group = new QVButtonGroup("Intensity scale:",box);
			intensity_mode_lin_ = new QRadioButton("Linear",intenstiy_group);
			intensity_mode_log_ = new QRadioButton("Log",intenstiy_group);
			grid->addWidget(box,1,1);	
			
			load();
		}
Spectrum3DCanvasPDP::~Spectrum3DCanvasPDP()
{
	
}
		
void Spectrum3DCanvasPDP::load()
{
	Spectrum3DCanvas* man = static_cast<Spectrum3DCanvas*>(manager_);
	
	if (man->getDotMode()==Spectrum3DCanvas::DOT_GRADIENT)
	{
		dot_mode_gradient_->setChecked(true);
	}
	else
	{
		if (man->getDotMode()==Spectrum3DCanvas::DOT_BLACK)
		{
			dot_mode_black_->setChecked(true);
		}
	}
	if (man->getShadeMode()==Spectrum3DCanvas::SHADE_FLAT)
	{
		shade_mode_flat_->setChecked(true);
	}
	else
	{
		if (man->getShadeMode()==Spectrum3DCanvas::SHADE_SMOOTH)
		{
			shade_mode_smooth_->setChecked(true);
		}
	}

	if(man->getIntScaleMode()==Spectrum3DCanvas::INT_LINEAR)
	{
		intensity_mode_lin_->setChecked(true);
	}
	else 
	{
		if (man->getIntScaleMode()==Spectrum3DCanvas::INT_LOG)
		{
			intensity_mode_log_->setChecked(true);
		}		
	}
	background_color_->setColor(QColor(man->getPrefAsString("Preferences:3D:BackgroundColor").c_str()));
	dot_gradient_->gradient().fromString(man->getPref("Preferences:3D:Dot:Gradient"));
	dot_interpolation_steps_->setValue(UnsignedInt(man->getPref("Preferences:3D:Dot:InterpolationSteps")));
	dot_line_width_->setValue(UnsignedInt(man->getPref("Preferences:3D:Dot:LineWidth")));
	axes_color_->setColor(QColor(man->getPrefAsString("Preferences:3D:AxesColor").c_str()));
}

void Spectrum3DCanvasPDP::save()
{
	Spectrum3DCanvas* man = static_cast<Spectrum3DCanvas*>(manager_);
	if(dot_mode_gradient_->isChecked())
	{
				man->setPref("Preferences:3D:Dot:Gradient",dot_gradient_->gradient().toString());
				man->setDotGradient(dot_gradient_->gradient().toString());

				man->setPref("Preferences:3D:Dot:InterpolationSteps",dot_interpolation_steps_->value());

				man->setPref("Preferences:3D:Dot:Mode",Spectrum3DCanvas::DOT_GRADIENT);
		
				if(shade_mode_flat_ -> isChecked())
				{
					man->setPref("Preferences:3D:Shade:Mode",Spectrum3DCanvas::SHADE_FLAT);
				}
				else if (shade_mode_smooth_->isChecked())
				{
					man->setPref("Preferences:3D:Shade:Mode",Spectrum3DCanvas::SHADE_SMOOTH);
				}
 			}	
	else
		{ 
			if(dot_mode_black_->isChecked())
				{
					man->setPref("Preferences:3D:Dot:Mode",Spectrum3DCanvas::DOT_BLACK);
				}
		} 
	if(intensity_mode_lin_->isChecked())
	{
		man->setPref("Preferences:3D:IntScale:Mode",Spectrum3DCanvas::INT_LINEAR);
	}
	else 
	{
		if(intensity_mode_log_->isChecked())
		{
			man->setPref("Preferences:3D:IntScale:Mode",Spectrum3DCanvas::INT_LOG);
			}
	}
	man->setPref("Preferences:3D:BackgroundColor",background_color_->getColor().name().ascii());
	man->setPref("Preferences:3D:AxesColor",axes_color_->getColor().name().ascii());
	man->setPref("Preferences:3D:Dot:LineWidth",dot_line_width_->value());
 	man->invalidate_();	
}
} // namespace Internal
} //namespace

