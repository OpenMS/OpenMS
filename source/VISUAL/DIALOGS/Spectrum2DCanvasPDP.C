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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DCanvasPDP.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/MultiGradientSelector.h>
#include <OpenMS/VISUAL/ColorSelector.h>

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
		
		Spectrum2DCanvasPDP::Spectrum2DCanvasPDP( Spectrum2DCanvas* manager, QWidget* parent, const char* name, WFlags f):
		PreferencesDialogPage(manager,parent,name,f)
		{
			help_ = "This is the preferences dialog of 2D spectrum!"
							"<br>";
		
			QGridLayout* grid;
			QLabel* label;
				
			grid = new QGridLayout(this,3,7);
			grid->setMargin(6);
			grid->setSpacing(4);	
			
			//dot mode
			QGroupBox* box = new QGroupBox(2,Qt::Horizontal,"Dot coloring",this);
			QVButtonGroup* coloring_group = new QVButtonGroup("Mode:",box);
			coloring_group->setFrameStyle(QFrame::NoFrame);
			box->addSpace(0);
			dot_mode_black_ = new QRadioButton("Black",coloring_group);
			dot_mode_gradient_ = new QRadioButton("Gradient",coloring_group);
			dot_gradient_ = new MultiGradientSelector(box);
			
			grid->addMultiCellWidget(box,0,1,0,0);
			
			//Surface mode
			box = new QGroupBox(2,Qt::Horizontal,"Surface coloring",this);
			surface_gradient_ = new MultiGradientSelector(box);
			grid->addWidget(box,0,1);

			//colors
			box = new QGroupBox(2,Qt::Horizontal,"Colors",this);
			label = new QLabel("Background color: ",box);
			background_color_ = new ColorSelector(box);
			label = new QLabel("Interpolation steps: ",box);
			interpolation_steps_ = new QSpinBox(10,1000,1,box,"");
			interpolation_steps_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			grid->addWidget(box,1,1);
			
			//details
			box = new QGroupBox(2,Qt::Horizontal,"Surface/contour details",this);
			label = new QLabel("Squares per axis: ",box);
			marching_squares_steps_ = new QSpinBox(10,100,1,box,"");
			marching_squares_steps_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			label = new QLabel("Contour lines: ",box);
			contour_steps_ = new QSpinBox(3,30,1,box,"");
			contour_steps_->setSizePolicy(QSizePolicy::Fixed,QSizePolicy::Minimum);
			grid->addWidget(box,2,0);

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
			man->setPref("Preferences:2D:BackgroundColor",background_color_->getColor().name().ascii());
			man->setPref("Preferences:2D:MarchingSquaresSteps",marching_squares_steps_->value());
			man->setPref("Preferences:2D:InterpolationSteps",interpolation_steps_->value());
			man->setPref("Preferences:2D:Contour:Lines",contour_steps_->value());
			
			man->repaintAll();
		}


	} // namespace Internal

} //namespace


