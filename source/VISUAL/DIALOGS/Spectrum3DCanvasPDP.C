// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/VISUAL/DIALOGS/Spectrum3DCanvasPDP.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>


// Qt
#include <QtGui/QLayout>
#include <QtGui/QRadioButton>
#include <QtGui/QLabel>
#include <QtGui/QGroupBox>
#include <QtGui/QSpinBox>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>
#include <QtGui/QButtonGroup>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
		Spectrum3DCanvasPDP::Spectrum3DCanvasPDP(Spectrum3DCanvas* manager, QWidget* parent)
				: PreferencesDialogPage(manager,parent)
		{
			help_ = "This is the preferences dialog of a 3D view of a map!";
			
			QGridLayout* grid = new QGridLayout(this);

			//peak color box
			QGroupBox* box = addBox(grid,0,0,"Peak colors",1,2);

			QVBoxLayout* tmp2 = new QVBoxLayout();
			dot_mode_black_= new QRadioButton("Black",this);
			tmp2->addWidget(dot_mode_black_);
			dot_mode_gradient_= new QRadioButton("Gradient",this);
			tmp2->addWidget(dot_mode_gradient_);
			addLayout(box->layout(),0,"Mode:",tmp2);

			dot_gradient_= new MultiGradientSelector(box);
			addWidget(box->layout(),1,"Gradient:",dot_gradient_);

			dot_interpolation_steps_= addSpinBox(box,10,1000,1);
			addWidget(box->layout(),2,"Interpolation steps:",dot_interpolation_steps_);
			finish(box->layout());
			
			tmp2 = new QVBoxLayout();
			QButtonGroup* tmp3 = new QButtonGroup(this);
			shade_mode_flat_= new QRadioButton("Flat",this);
			tmp2->addWidget(shade_mode_flat_);
			tmp3->addButton(shade_mode_flat_);
			shade_mode_smooth_= new QRadioButton("Smooth",this);
			tmp2->addWidget(shade_mode_smooth_);
			tmp3->addButton(shade_mode_smooth_);
			addLayout(box->layout(),3,"Shade mode:",tmp2);
			finish(box->layout());
			
			//misc box
			box = addBox(grid,1,0,"Misc");
			
			dot_line_width_ = addSpinBox(box,1,10,1);
			addWidget(box->layout(),0,"Line width:",dot_line_width_);
			
			background_color_= new ColorSelector(box);
			addWidget(box->layout(),1,"Background color:",background_color_);
			axes_color_= new ColorSelector(box);
			addWidget(box->layout(),2,"Axis color:",axes_color_);
			finish(box->layout());

			//data reduction box
		 	box = addBox(grid,1,1,"Data reduction");
			data_reduction_= new QComboBox( box);
			data_reduction_->insertItem(0,"Reduction OFF");
			data_reduction_->insertItem(1,"MaxReduction");
			data_reduction_->insertItem(2,"SumReduction");
			addWidget(box->layout(),0,"Mode:",data_reduction_);

			reduction_diplay_peaks_= addSpinBox(box,5000,200000,5000);
			addWidget(box->layout(),1,"Displayed Peaks:",reduction_diplay_peaks_);
			finish(box->layout());
			
			finish(grid);
		
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
			data_reduction_->setCurrentIndex(data_reduction_->findText(man->getPrefAsString("Preferences:3D:Reduction:Mode").c_str()));	
			reduction_diplay_peaks_->setValue(UnsignedInt(man->getPrefAsInt("Preferences:3D:DisplayedPeaks")));
		
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

			if(data_reduction_->currentText().toAscii().data()=="MaxReduction")
			{
				man->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_MAX);
			}
			else if(data_reduction_->currentText().toAscii().data()=="Reduction OFF")
			{
				man->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_OFF);
			}
			else if(data_reduction_->currentText().toAscii().data()=="SumReduction")
			{
				man->setPref("Preferences:3D:Data:Mode",Spectrum3DCanvas::REDUCTION_SUM);
			}
			
			man->setPref("Preferences:3D:Reduction:Mode", data_reduction_->currentText().toAscii().data());
			man->setPref("Preferences:3D:DisplayedPeaks",	reduction_diplay_peaks_->value());
			
			man->setPref("Preferences:3D:BackgroundColor",background_color_->getColor().name().toAscii().data());
			man->setPref("Preferences:3D:AxesColor",axes_color_->getColor().name().toAscii().data());
			man->setPref("Preferences:3D:Dot:LineWidth",dot_line_width_->value());
		
			man->setDataMode();
			man->repaintAll();	
		

		}
	} // namespace Internal
} //namespace

