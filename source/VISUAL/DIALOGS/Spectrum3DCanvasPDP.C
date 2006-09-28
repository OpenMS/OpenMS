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
#include <qhbuttongroup.h>
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
			grid = new QGridLayout(this, 4, 2);
			grid->setMargin(6);
			grid->setSpacing(4);	

			QGroupBox* box = new QGroupBox(2,Qt::Horizontal,"Dot coloring",this);

			QVButtonGroup* coloring_group = new QVButtonGroup("Color Mode:",box);
			box->addSpace(0);  
			dot_mode_black_ = new QRadioButton("Black",coloring_group);
			dot_mode_gradient_ = new QRadioButton("Gradient",coloring_group);
			dot_gradient_ = new MultiGradientSelector(coloring_group);
		// 	box->addSpace(0);
			//		
			QHButtonGroup* interpolation_box = new QHButtonGroup("Interpolation steps",box);
			label = new QLabel("Interpolation steps: ",interpolation_box);
			dot_interpolation_steps_ = new QSpinBox(10,1000,1,interpolation_box,"");
			dot_interpolation_steps_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			box->addSpace(0);
			QVButtonGroup* shading_group = new QVButtonGroup("Shade Mode:",box);
			//shading_group->setFrameStyle(QFrame::NoFrame);
			shade_mode_flat_ = new QRadioButton("Flat",shading_group);
			shade_mode_smooth_ = new QRadioButton("Smooth",shading_group);	
			grid->addMultiCellWidget(box,0,3,0,0);

			box = new QGroupBox(2,Qt::Horizontal,"Line Width",this);
			label = new QLabel("Line Width: ",box);
			dot_line_width_ = new QSpinBox(1,10,1,box,"");
			dot_line_width_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			grid->addWidget(box,0,1);	

			box = new QGroupBox(2,Qt::Horizontal,"Colors",this);
			label = new QLabel("Background color: ",box);
			background_color_ = new ColorSelector(box);
			label = new QLabel("Axes Color: ",box);
			axes_color_ = new ColorSelector(box);
			grid->addWidget(box,1,1);	

			box = new QGroupBox(2,Qt::Horizontal,"Data",this);
			QVButtonGroup* data_group = new QVButtonGroup("DataReduction",box);
			reduction_off_ = new QRadioButton("Off",data_group);	
			reduction_on_max_ = new QRadioButton("Max-Reduction",data_group);
			reduction_on_sum_ = new QRadioButton("Sum-Reduction",data_group);			
			box->addSpace(0);
			
			label = new QLabel("Number of peaks per Reductionstep:(Max-Red.) ",box);
			reduction_ratio_max_ = new QSpinBox(10,100,1,box,"");
		
			label = new QLabel("m/z-Range per Reductionstep:(Sum-Red.) ",box);
			reduction_ratio_sum_ = new QSpinBox(10,100,1,box,"");
			grid->addWidget(box,2,1);	

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

			if (man->getDataMode()==Spectrum3DCanvas::REDUCTION_MAX)
			{
				reduction_on_max_->setChecked(true);
			}
			else if (man->getDataMode()==Spectrum3DCanvas::REDUCTION_OFF)
			{
				reduction_off_->setChecked(true);
			}
			else if(man->getDataMode()==Spectrum3DCanvas::REDUCTION_SUM)
			{
				reduction_on_sum_->setChecked(true);
				
			}
			reduction_ratio_max_->setValue(UnsignedInt(manager_->getPref("Preferences:3D:Data:Reduction:Max")));
			reduction_ratio_sum_->setValue(UnsignedInt(manager_->getPref("Preferences:3D:Data:Reduction:Sum")));

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
			if(reduction_on_max_->isChecked())
			{
				man->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_MAX);
			}
			else if(reduction_off_->isChecked())
			{
				man->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_OFF);
			}
			else if(reduction_on_sum_->isChecked())
			{
				man->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_SUM);
			}
			man->setPref("Preferences:3D:Data:Reduction:Max",reduction_ratio_max_->value());
			man->setPref("Preferences:3D:Data:Reduction:Sum",reduction_ratio_sum_->value());
			man->setDataMode(man->getPrefAsInt("Preferences:3D:Data:Mode"));
			
			man->setPref("Preferences:3D:BackgroundColor",background_color_->getColor().name().ascii());
			man->setPref("Preferences:3D:AxesColor",axes_color_->getColor().name().ascii());
			man->setPref("Preferences:3D:Dot:LineWidth",dot_line_width_->value());
			
			man->repaintAll();
		}
	} // namespace Internal
} //namespace

