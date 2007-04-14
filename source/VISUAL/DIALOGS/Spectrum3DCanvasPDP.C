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
			grid->setMargin(0);
			
			//peak color box
			QGroupBox* box =	addBox_(grid,0,0,"Peak colors",1,2);

			dot_gradient_= new MultiGradientSelector(box);
			addWidget_(box->layout(),0,"Gradient:",dot_gradient_);
			finish_(box->layout());
			
			QVBoxLayout* tmp2 = new QVBoxLayout();
			QButtonGroup* tmp3 = new QButtonGroup(this);
			shade_mode_flat_= new QRadioButton("Flat",this);
			tmp2->addWidget(shade_mode_flat_);
			tmp3->addButton(shade_mode_flat_);
			shade_mode_smooth_= new QRadioButton("Smooth",this);
			tmp2->addWidget(shade_mode_smooth_);
			tmp3->addButton(shade_mode_smooth_);
			addLayout_(box->layout(),2,"Shade mode:",tmp2);
			finish_(box->layout());
			
			//misc box
			box =	addBox_(grid,1,0,"Misc");
			
			dot_line_width_ = addSpinBox_(box,1,10,1);
			addWidget_(box->layout(),0,"Line width:",dot_line_width_);
			
			finish_(box->layout());

			//data reduction box
		 	box =	addBox_(grid,1,1,"Data reduction");
			data_reduction_= new QComboBox( box);
			data_reduction_->insertItem(0,"Reduction OFF");
			data_reduction_->insertItem(1,"MaxReduction");
			data_reduction_->insertItem(2,"SumReduction");
			addWidget_(box->layout(),0,"Mode:",data_reduction_);

			reduction_diplay_peaks_= addSpinBox_(box,5000,200000,5000);
			addWidget_(box->layout(),1,"Displayed Peaks:",reduction_diplay_peaks_);
			finish_(box->layout());
			
			finish_(grid);
		
			load();
		}

		Spectrum3DCanvasPDP::~Spectrum3DCanvasPDP()
		{
			
		}
				
		void Spectrum3DCanvasPDP::load()
		{
			Spectrum3DCanvas* man = static_cast<Spectrum3DCanvas*>(manager_);
			
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
			reduction_diplay_peaks_->setValue(UInt(man->getPrefAsInt("Preferences:3D:DisplayedPeaks")));

			dot_gradient_->gradient().fromString(man->getPref("Preferences:3D:Dot:Gradient"));
			dot_line_width_->setValue(UInt(man->getPref("Preferences:3D:Dot:LineWidth")));
		}
		
		void Spectrum3DCanvasPDP::save()
		{
			Spectrum3DCanvas* man = static_cast<Spectrum3DCanvas*>(manager_);
			man->setPref("Preferences:3D:Dot:Gradient",dot_gradient_->gradient().toString());
			man->setDotGradient(dot_gradient_->gradient().toString());
	
			if(shade_mode_flat_ -> isChecked())
			{
				man->setPref("Preferences:3D:Shade:Mode",Spectrum3DCanvas::SHADE_FLAT);
			}
			else if (shade_mode_smooth_->isChecked())
			{
				man->setPref("Preferences:3D:Shade:Mode",Spectrum3DCanvas::SHADE_SMOOTH);
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
			
			man->setPref("Preferences:3D:Dot:LineWidth",dot_line_width_->value());
		
			man->repaintAll();	
		

		}
	} // namespace Internal
} //namespace

