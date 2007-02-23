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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DCanvasPDP.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/ColorSelector.h>

// Qt
#include <QtGui/QLayout>
#include <QtGui/QRadioButton>
#include <QtGui/QLabel>
#include <QtGui/QGroupBox>
#include <QtGui/QSpinBox>
#include <QtGui/QGridLayout>
#include <QtGui/QButtonGroup>

using namespace std;

namespace OpenMS
{

	namespace Internal
	{
		
		Spectrum2DCanvasPDP::Spectrum2DCanvasPDP( Spectrum2DCanvas* manager, QWidget* parent):
		PreferencesDialogPage(manager,parent)
		{
			help_ = "This is the preferences dialog of a 2D view of a map!";
		
			QGridLayout* grid = new QGridLayout(this);
			
			//colors
			QGroupBox* box = addBox(grid,0,0,"Colors",2,1);
			background_color_ = new ColorSelector(box);
			addWidget(box->layout(),0,"Background color:",background_color_);
			interpolation_steps_ = addSpinBox(box,10,1000,1);
			addWidget(box->layout(),1,"Interpolation steps:",interpolation_steps_);
			finish(box->layout());		
						
			//dot mode
			box = addBox(grid,0,1,"Dot Colors");
			QVBoxLayout* tmp2 = new QVBoxLayout();
			dot_mode_black_ = new QRadioButton("Black",this);
			tmp2->addWidget(dot_mode_black_);
			dot_mode_gradient_ = new QRadioButton("Gradient",this);
			tmp2->addWidget(dot_mode_gradient_);
			addLayout(box->layout(),0,"Mode:",tmp2);
	
			dot_gradient_ = new MultiGradientSelector(box);
			addWidget(box->layout(),1,"Gradient:",dot_gradient_);
			finish(box->layout());
					
			//surface mode
			box = addBox(grid,1,1,"Surface/contour settings");
			surface_gradient_ = new MultiGradientSelector(box);
			addWidget(box->layout(),0,"Gradient:",surface_gradient_);			

			marching_squares_steps_ = addSpinBox(box,10,100,1);
			addWidget(box->layout(),1,"Squares per axis:",marching_squares_steps_);
			
			contour_steps_ = addSpinBox(box,3,30,1);
			addWidget(box->layout(),2,"Contour lines:",contour_steps_);
			finish(box->layout());

			load();
		}
					
		Spectrum2DCanvasPDP::~Spectrum2DCanvasPDP()
		{
			
		}
		
		void Spectrum2DCanvasPDP::load()
		{
			Spectrum2DCanvas* man = static_cast<Spectrum2DCanvas*>(manager_);
			
			if (man->getDotMode()==Spectrum2DCanvas::DOT_GRADIENT)
			{
				dot_mode_gradient_->setChecked(true);
			}
			else if (man->getDotMode()==Spectrum2DCanvas::DOT_BLACK)
			{
				dot_mode_black_->setChecked(true);
			}
			dot_gradient_->gradient().fromString(man->getPrefAsString("Preferences:2D:Dot:Gradient"));
			surface_gradient_->gradient().fromString(man->getPrefAsString("Preferences:2D:Surface:Gradient"));
			background_color_->setColor(QColor(man->getPrefAsString("Preferences:2D:BackgroundColor").c_str()));
			marching_squares_steps_->setValue(UnsignedInt(man->getPref("Preferences:2D:MarchingSquaresSteps")));
			interpolation_steps_->setValue(UnsignedInt(man->getPref("Preferences:2D:InterpolationSteps")));
			contour_steps_->setValue(UnsignedInt(man->getPref("Preferences:2D:Contour:Lines")));
		}
		
		void Spectrum2DCanvasPDP::save()
		{
			Spectrum2DCanvas* man = static_cast<Spectrum2DCanvas*>(manager_);
			
			if (dot_mode_gradient_->isChecked())
			{
				man->setDotMode(Spectrum2DCanvas::DOT_GRADIENT);
			}
			else if (dot_mode_black_->isChecked())
			{
				man->setDotMode(Spectrum2DCanvas::DOT_BLACK);
			}
			man->setDotGradient(dot_gradient_->gradient().toString());
			man->setSurfaceGradient(surface_gradient_->gradient().toString());
			man->setPref("Preferences:2D:BackgroundColor",background_color_->getColor().name().toAscii().data());
			man->setPref("Preferences:2D:MarchingSquaresSteps",marching_squares_steps_->value());
			man->setPref("Preferences:2D:InterpolationSteps",interpolation_steps_->value());
			man->setPref("Preferences:2D:Contour:Lines",contour_steps_->value());
			
			man->repaintAll();
		}


	} // namespace Internal

} //namespace


